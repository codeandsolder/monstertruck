//! Physical line and circular-arc coverage, independent of edge subdivision.
//! Samples propose a circle; whole-domain rational certificates prove it.
use super::{
    ANALYTIC_BOUNDARY_PROOF_TOL_MM, BSplineSupport, BoundaryChainKey, CellKind, CellWitnesses,
    CircleLocusKey, CurveSupport, DIRECTION_Q, EdgeGeometry, OrientedSupportKey, POINT_Q_MM,
    PointKey, SupportLocusKey, brep, bspline_domain, bspline_line_support,
    canonical_boundary_direction, circle_lies_on_support, circle_section_of_plane_and_support,
    dequantize_scalar, dequantize_vector, evaluate_bspline, line_lies_on_support,
    line_segment_lies_on_cone, support_locus_key, support_pair_can_intersect,
};
use crate::math3::{add, cross, distance, dot, mul, norm, normalize, sub};
use monstertruck_geometry::{
    analytic_carrier::AnalyticSupport,
    analytic_curve_incidence::{
        CircularArcFrame, CurveIncidence, certify_circular_arc_coverage,
        certify_cylindrical_arc_coverage, certify_nurbs_curve_on_support,
    },
    analytic_line_segment::{SegmentLocus, certify_nurbs_segment_locus},
    prelude::{BoundedCurve, NurbsCurve, ParametricCurve, Point3, Vector3, Vector4},
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

const TOL: f64 = ANALYTIC_BOUNDARY_PROOF_TOL_MM;

const fn point(p: [f64; 3]) -> Point3 { Point3::new(p[0], p[1], p[2]) }
const fn vector(v: [f64; 3]) -> Vector3 { Vector3::new(v[0], v[1], v[2]) }
fn decode(p: PointKey) -> Option<[f64; 3]> { dequantize_vector(p, POINT_Q_MM) }

fn sphere_or_torus(support: &SupportLocusKey) -> Option<AnalyticSupport> {
    match support {
        SupportLocusKey::Sphere {
            center_mm,
            radius_mm,
        } => Some(AnalyticSupport::Sphere {
            center: point(decode(*center_mm)?),
            radius: dequantize_scalar(*radius_mm)?,
        }),
        SupportLocusKey::Torus {
            center_mm,
            axis,
            major_radius_mm,
            minor_radius_mm,
        } => Some(AnalyticSupport::Torus {
            center: point(decode(*center_mm)?),
            axis: vector(normalize(dequantize_vector(*axis, DIRECTION_Q)?, 0.)?),
            major_radius: dequantize_scalar(*major_radius_mm)?,
            minor_radius: dequantize_scalar(*minor_radius_mm)?,
        }),
        _ => None,
    }
}

/// Four exact rational quadrants cover the whole untrimmed circle.
fn full_circle(circle: brep::CircleSupport) -> Option<NurbsCurve<Vector4>> {
    let normal = normalize(circle.normal, 0.)?;
    let x = normalize(
        sub(
            circle.x_direction,
            mul(normal, dot(normal, circle.x_direction)),
        ),
        0.,
    )?;
    let y = cross(normal, x);
    let rays = [x, y, mul(x, -1.), mul(y, -1.), x];
    let mut controls = Vec::with_capacity(9);
    for pair in rays.windows(2) {
        controls.push(add(circle.center_mm, mul(pair[0], circle.radius_mm)));
        controls.push(add(
            circle.center_mm,
            mul(add(pair[0], pair[1]), circle.radius_mm),
        ));
    }
    controls.push(add(circle.center_mm, mul(x, circle.radius_mm)));
    BSplineSupport {
        degree: 2,
        control_points_mm: controls,
        knots: vec![0., 0., 0., 1., 1., 2., 2., 3., 3., 4., 4., 4.],
        weights: Some(
            (0..9)
                .map(|i| {
                    if i % 2 == 0 {
                        1.
                    } else {
                        std::f64::consts::FRAC_1_SQRT_2
                    }
                })
                .collect(),
        ),
    }
    .rational_curve()
}

/// Exact range of distance from every point of a circle to a sphere center.
/// A nonconcentric circular section is allowed; endpoint samples alone
/// cannot determine the radial extremum.
pub(super) fn circle_on_sphere(circle: brep::CircleSupport, center: [f64; 3], radius: f64) -> bool {
    let delta = sub(circle.center_mm, center);
    let height = dot(delta, circle.normal);
    let radial = norm(sub(delta, mul(circle.normal, height)));
    let inner = height.hypot(circle.radius_mm - radial);
    let outer = height.hypot(circle.radius_mm + radial);
    (inner - radius).abs().max((outer - radius).abs()) <= TOL
}

/// Rotating a circle chart through the acute angle between normals moves
/// every point by at most sqrt(2)*radius*sin(angle). This bound follows
/// from 2*sin(angle/2) <= sqrt(2)*sin(angle) for 0 <= angle <= pi/2.
fn plane_rotation_bound(radius: f64, first: [f64; 3], second: [f64; 3]) -> f64 {
    std::f64::consts::SQRT_2 * radius * norm(cross(first, second))
}

/// Compare against an exact coaxial parallel or meridian of the ring torus.
/// The distance-to-tube-center-circle function is 1-Lipschitz; translating
/// the center and rotating the chart have additive global error bounds.
/// Unlike the conservative quartic coefficient test, this handles small
/// quantization errors at tangent parallels without losing tightness.
fn circle_on_torus(
    circle: brep::CircleSupport,
    center: [f64; 3],
    axis: [f64; 3],
    major: f64,
    minor: f64,
) -> bool {
    if !(major > minor && minor > TOL) {
        return false;
    }
    let delta = sub(circle.center_mm, center);
    let height = dot(delta, axis);
    let radial = sub(delta, mul(axis, height));
    let radial_length = norm(radial);
    let parallel_error = radial_length
        + plane_rotation_bound(circle.radius_mm, circle.normal, axis)
        + ((circle.radius_mm - major).hypot(height) - minor).abs();
    if parallel_error <= TOL {
        return true;
    }
    normalize(radial, 0.).is_some_and(|direction| {
        let meridian_normal = cross(axis, direction);
        let meridian_error = (radial_length - major).hypot(height)
            + plane_rotation_bound(circle.radius_mm, circle.normal, meridian_normal)
            + (circle.radius_mm - minor).abs();
        meridian_error <= TOL
    })
}

pub(super) fn circle_on_curved_support(
    circle: brep::CircleSupport,
    support: &SupportLocusKey,
) -> bool {
    let Some(normal) = normalize(circle.normal, 0.) else {
        return false;
    };
    if !circle.radius_mm.is_finite() || circle.radius_mm <= TOL {
        return false;
    }
    let circle = brep::CircleSupport { normal, ..circle };
    match support {
        SupportLocusKey::Sphere {
            center_mm,
            radius_mm,
        } => decode(*center_mm)
            .zip(dequantize_scalar(*radius_mm))
            .is_some_and(|(center, radius)| circle_on_sphere(circle, center, radius)),
        SupportLocusKey::Torus {
            center_mm,
            axis,
            major_radius_mm,
            minor_radius_mm,
        } => {
            let torus = || {
                Some((
                    decode(*center_mm)?,
                    normalize(dequantize_vector(*axis, DIRECTION_Q)?, 0.)?,
                    dequantize_scalar(*major_radius_mm)?,
                    dequantize_scalar(*minor_radius_mm)?,
                ))
            };
            if torus().is_some_and(|(center, axis, major, minor)| {
                circle_on_torus(circle, center, axis, major, minor)
            }) {
                return true;
            }
            sphere_or_torus(support)
                .zip(full_circle(circle))
                .is_some_and(|(support, curve)| {
                    certify_nurbs_curve_on_support(&curve, support, TOL) == CurveIncidence::Proven
                })
        }
        _ => false,
    }
}

fn on_circle(p: [f64; 3], circle: brep::CircleSupport) -> bool {
    let delta = sub(p, circle.center_mm);
    dot(delta, circle.normal).abs() <= TOL && (norm(delta) - circle.radius_mm).abs() <= TOL
}
pub(super) fn on_line(p: [f64; 3], line: brep::LineSupport) -> bool {
    normalize(line.direction, 0.)
        .is_some_and(|direction| norm(cross(sub(p, line.origin_mm), direction)) <= TOL)
}
fn line_incidence(
    edge: &EdgeGeometry,
    line: brep::LineSupport,
    loci: &[SupportLocusKey; 2],
) -> bool {
    edge.start != edge.end
        && [edge.start, edge.end]
            .into_iter()
            .all(|p| decode(p).is_some_and(|p| on_line(p, line)))
        && loci.iter().all(|locus| match locus {
            SupportLocusKey::Cone { .. } => {
                line_segment_lies_on_cone(line, [edge.start, edge.end], locus)
            }
            _ => line_lies_on_support(line, locus),
        })
}

fn fit_circle(spline: &BSplineSupport) -> Option<brep::CircleSupport> {
    let (lo, hi) = bspline_domain(spline)?;
    let first = evaluate_bspline(spline, lo)?;
    let second = evaluate_bspline(spline, (hi - lo).mul_add(1. / 3., lo))?;
    let third = evaluate_bspline(spline, (hi - lo).mul_add(2. / 3., lo))?;
    let u = sub(second, first);
    let v = sub(third, first);
    let normal = cross(u, v);
    let denominator = 2. * dot(normal, normal);
    if !denominator.is_finite() || denominator <= 0. {
        return None;
    }
    let displacement = mul(
        add(
            mul(cross(v, normal), dot(u, u)),
            mul(cross(normal, u), dot(v, v)),
        ),
        1. / denominator,
    );
    let center = add(first, displacement);
    let radius = norm(displacement);
    if !radius.is_finite() || radius <= TOL {
        return None;
    }
    Some(brep::CircleSupport {
        center_mm: center,
        radius_mm: radius,
        normal: normalize(normal, 0.)?,
        x_direction: normalize(sub(first, center), 0.)?,
    })
}

fn positive_angle(start: [f64; 3], end: [f64; 3], normal: [f64; 3]) -> f64 {
    dot(normal, cross(start, end))
        .atan2(dot(start, end))
        .rem_euclid(TAU)
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Arc {
    start: PointKey,
    end: PointKey,
    sweep: f64,
}

pub(super) enum BoundedEdge {
    Line,
    Circle(CircleLocusKey, Arc),
    CylinderSection(Arc),
    Spline(super::spline_strip::ProfilePiece),
}

fn circle_arc(
    edge: &EdgeGeometry,
    circle: brep::CircleSupport,
    spline: Option<&NurbsCurve<Vector4>>,
    loci: &[SupportLocusKey; 2],
    cells: &mut CellWitnesses,
) -> Option<BoundedEdge> {
    if !circle.radius_mm.is_finite() || circle.radius_mm <= TOL {
        return None;
    }
    let raw_normal = normalize(circle.normal, 0.)?;
    let (normal, normal_key) = canonical_boundary_direction(raw_normal)?;
    let circle = brep::CircleSupport { normal, ..circle };
    if !loci
        .iter()
        .all(|locus| circle_lies_on_support(circle, locus))
    {
        return None;
    }
    let start = decode(edge.start)?;
    let end = decode(edge.end)?;
    if !on_circle(start, circle) || !on_circle(end, circle) {
        return None;
    }
    let tol_rad = TOL / circle.radius_mm;
    let mut ends = [edge.start, edge.end];
    let sweep = if let Some(spline) = spline {
        let frame = CircularArcFrame::new(
            point(circle.center_mm),
            vector(normal),
            circle.radius_mm,
            point(start),
            TOL,
        )?;
        let coverage = certify_circular_arc_coverage(spline, frame)?;
        let near = |a: Point3, b: [f64; 3]| distance([a.x, a.y, a.z], b) <= TOL;
        let from = frame.point(coverage.start_rad);
        let to = frame.point(coverage.start_rad + coverage.sweep_rad);
        if near(from, start) && near(to, end) {
            // Already increasing in canonical direction.
        } else if near(from, end) && near(to, start) {
            ends.reverse();
        } else {
            return None;
        }
        coverage.sweep_rad
    } else {
        if edge.parameter_forward != (dot(raw_normal, normal) > 0.) {
            ends.reverse();
        }
        if edge.start == edge.end {
            if edge.start_vertex != edge.end_vertex {
                return None;
            }
            TAU
        } else {
            positive_angle(
                sub(decode(ends[0])?, circle.center_mm),
                sub(decode(ends[1])?, circle.center_mm),
                normal,
            )
        }
    };
    if !sweep.is_finite()
        || sweep <= tol_rad
        || sweep > TAU + tol_rad
        || (ends[0] == ends[1]) != ((sweep - TAU).abs() <= tol_rad)
    {
        return None;
    }
    cells.record_direction(normal_key, normal);
    let locus = CircleLocusKey {
        center_mm: cells.point(circle.center_mm)?,
        radius_mm: cells.scalar(CellKind::Radius, circle.radius_mm)?,
        normal: normal_key,
    };
    Some(BoundedEdge::Circle(
        locus,
        Arc {
            start: ends[0],
            end: ends[1],
            sweep,
        },
    ))
}

/// The axial projection of a nonparallel plane/cylinder intersection
/// is bijective onto its circle. Establish full physical plane incidence
/// and full cylindrical/angular coverage independently. Bernstein control
/// rays are measured in the transverse chart of the actual rational curve,
/// without sampling, refitting or allocating a projected curve.
struct CylinderSectionChart {
    origin: [f64; 3],
    axis: [f64; 3],
    radius: f64,
    plane_origin: [f64; 3],
    plane_normal: [f64; 3],
}

impl CylinderSectionChart {
    fn new(loci: &[SupportLocusKey; 2]) -> Option<Self> {
        let (normal, offset, axis, closest, radius) = match (&loci[0], &loci[1]) {
            (
                SupportLocusKey::Plane { normal, offset_mm },
                SupportLocusKey::Cylinder {
                    axis,
                    closest_mm,
                    radius_mm,
                },
            )
            | (
                SupportLocusKey::Cylinder {
                    axis,
                    closest_mm,
                    radius_mm,
                },
                SupportLocusKey::Plane { normal, offset_mm },
            ) => (*normal, *offset_mm, *axis, *closest_mm, *radius_mm),
            _ => return None,
        };
        let plane = dequantize_vector(normal, DIRECTION_Q)?;
        let plane_normal = normalize(plane, 0.)?;
        let axis = normalize(dequantize_vector(axis, DIRECTION_Q)?, 0.)?;
        // Transverse sections already have a true circle certificate;
        // axial planes can intersect the cylinder in one or two lines.
        if dot(plane_normal, axis).abs() <= 4. * DIRECTION_Q
            || norm(cross(plane_normal, axis)) <= 4. * DIRECTION_Q
        {
            return None;
        }
        let radius = dequantize_scalar(radius)?;
        if radius <= TOL {
            return None;
        }
        Some(Self {
            origin: decode(closest)?,
            axis,
            radius,
            plane_origin: mul(plane, dequantize_scalar(offset)? / dot(plane, plane)),
            plane_normal,
        })
    }

    fn project(&self, p: [f64; 3]) -> [f64; 3] {
        let delta = sub(p, self.origin);
        sub(delta, mul(self.axis, dot(delta, self.axis)))
    }

    fn on_section(&self, p: [f64; 3]) -> bool {
        dot(sub(p, self.plane_origin), self.plane_normal).abs() <= TOL
            && (norm(self.project(p)) - self.radius).abs() <= TOL
    }

    fn certify(&self, edge: &EdgeGeometry, curve: &NurbsCurve<Vector4>) -> Option<Arc> {
        let plane = AnalyticSupport::Plane {
            origin: point(self.plane_origin),
            normal: vector(self.plane_normal),
        };
        if certify_nurbs_curve_on_support(curve, plane, TOL) != CurveIncidence::Proven {
            return None;
        }
        let start = decode(edge.start)?;
        let end = decode(edge.end)?;
        if !self.on_section(start) || !self.on_section(end) {
            return None;
        }
        let (lo, hi) = curve.range_tuple();
        let near = |p: Point3, q: [f64; 3]| distance([p.x, p.y, p.z], q) <= TOL;
        let a = curve.evaluate(lo);
        let b = curve.evaluate(hi);
        if !((near(a, start) && near(b, end)) || (near(a, end) && near(b, start))) {
            return None;
        }
        let projected_start = add(self.origin, self.project(start));
        let projected_end = add(self.origin, self.project(end));
        let frame = CircularArcFrame::new(
            point(self.origin),
            vector(self.axis),
            self.radius,
            point(projected_start),
            TOL,
        )?;
        let cover = certify_cylindrical_arc_coverage(curve, frame)?;
        let from = frame.point(cover.start_rad);
        let to = frame.point(cover.start_rad + cover.sweep_rad);
        let ends = if near(from, projected_start) && near(to, projected_end) {
            [edge.start, edge.end]
        } else if near(from, projected_end) && near(to, projected_start) {
            [edge.end, edge.start]
        } else {
            return None;
        };
        let tol = TOL / self.radius;
        let sweep = cover.sweep_rad;
        if !sweep.is_finite()
            || sweep <= tol
            || sweep > TAU + tol
            || (ends[0] == ends[1]) != ((sweep - TAU).abs() <= tol)
            || (ends[0] == ends[1] && edge.start_vertex != edge.end_vertex)
        {
            return None;
        }
        Some(Arc {
            start: ends[0],
            end: ends[1],
            sweep,
        })
    }
}

/// A declared nonsimple curve may still be proved to traverse this
/// ellipse exactly once. No general curve recognition relies on this path.
pub(super) fn certify_flagged_cylinder_section(
    edge: &EdgeGeometry,
    curve: &CurveSupport,
    loci: &[SupportLocusKey; 2],
) -> Option<BoundedEdge> {
    let CurveSupport::BSpline(spline) = curve else {
        return None;
    };
    let chart = CylinderSectionChart::new(loci)?;
    chart
        .certify(edge, &spline.rational_curve()?)
        .map(BoundedEdge::CylinderSection)
}

pub(super) fn classify(
    edge: &EdgeGeometry,
    curve: &CurveSupport,
    loci: &[SupportLocusKey; 2],
    cells: &mut CellWitnesses,
) -> Option<BoundedEdge> {
    if !support_pair_can_intersect(loci) {
        return None;
    }
    match curve {
        CurveSupport::Line(line) => line_incidence(edge, *line, loci).then_some(BoundedEdge::Line),
        CurveSupport::Circle(circle) => circle_arc(edge, *circle, None, loci, cells),
        CurveSupport::BSpline(spline) => {
            let rational = spline.rational_curve()?;
            if let Some(line) = bspline_line_support(spline)
                && line_incidence(edge, line, loci)
                && certify_nurbs_segment_locus(
                    &rational,
                    point(decode(edge.start)?),
                    point(decode(edge.end)?),
                    TOL,
                ) == SegmentLocus::WithinTolerance
            {
                return Some(BoundedEdge::Line);
            }
            if let Some(arc) =
                CylinderSectionChart::new(loci).and_then(|chart| chart.certify(edge, &rational))
            {
                return Some(BoundedEdge::CylinderSection(arc));
            }
            if let Some(key) = super::spline_strip::certified_edge(edge, spline, loci) {
                return Some(BoundedEdge::Spline(key));
            }
            let proposed = circle_section_of_plane_and_support(&loci[0], &loci[1])
                .or_else(|| circle_section_of_plane_and_support(&loci[1], &loci[0]))
                .and_then(|(center, normal, radius)| {
                    Some(brep::CircleSupport {
                        center_mm: center,
                        normal,
                        radius_mm: radius,
                        x_direction: normalize(sub(decode(edge.start)?, center), 0.)?,
                    })
                });
            proposed
                .and_then(|circle| circle_arc(edge, circle, Some(&rational), loci, cells))
                .or_else(|| circle_arc(edge, fit_circle(spline)?, Some(&rational), loci, cells))
        }
        CurveSupport::Other { .. } => None,
    }
}

/// Directed chains retain arc choice, unlike sorted endpoint pairs. Reject
/// branch, repeated coverage and winding beyond one full revolution.
pub(super) fn normalize_arcs(
    supports: &[OrientedSupportKey; 2],
    locus: &CircleLocusKey,
    arcs: &[Arc],
) -> Result<Vec<BoundaryChainKey>, String> {
    let radius = dequantize_scalar(locus.radius_mm).ok_or("invalid circle radius")?;
    let tol = TOL / radius;
    let center = decode(locus.center_mm).ok_or("invalid circle center")?;
    let normal = normalize(
        dequantize_vector(locus.normal, DIRECTION_Q).ok_or("invalid circle normal")?,
        0.,
    )
    .ok_or("invalid circle direction")?;
    normalize_angular_chains(arcs, tol, normal, |p| sub(p, center)).map(|chains| {
        chains
            .into_iter()
            .map(|endpoints_mm| BoundaryChainKey::Circular {
                supports: supports.clone(),
                locus: locus.clone(),
                endpoints_mm,
            })
            .collect()
    })
}

pub(super) fn normalize_cylinder_sections(
    supports: &[OrientedSupportKey; 2],
    arcs: &[Arc],
) -> Result<Vec<BoundaryChainKey>, String> {
    let chart = CylinderSectionChart::new(&supports.each_ref().map(support_locus_key))
        .ok_or("invalid oblique cylinder section")?;
    normalize_angular_chains(arcs, TOL / chart.radius, chart.axis, |p| chart.project(p)).map(
        |chains| {
            chains
                .into_iter()
                .map(|endpoints_mm| BoundaryChainKey::CylinderSection {
                    supports: supports.clone(),
                    endpoints_mm,
                })
                .collect()
        },
    )
}

/// Chain topology stays in physical 3-D coordinates; only angular coverage
/// is evaluated in a transverse chart. This preserves actual boundary
/// vertices and rejects repeated coverage regardless of curve subdivision.
fn normalize_angular_chains(
    arcs: &[Arc],
    tol: f64,
    normal: [f64; 3],
    radial: impl Fn([f64; 3]) -> [f64; 3],
) -> Result<Vec<Option<[PointKey; 2]>>, String> {
    let mut outgoing = BTreeMap::new();
    let mut incoming = BTreeSet::new();
    for (i, arc) in arcs.iter().enumerate() {
        if outgoing.insert(arc.start, i).is_some() || !incoming.insert(arc.end) {
            return Err("circular boundary branches or overlaps at endpoint".into());
        }
    }
    // Open chains first, so partition vertices never become endpoints.
    let mut seeds = (0..arcs.len()).collect::<Vec<_>>();
    seeds.sort_by_key(|&i| incoming.contains(&arcs[i].start));
    let mut visited = BTreeSet::new();
    let mut result = Vec::new();
    let mut intervals = Vec::new();
    let angle_origin = arcs
        .first()
        .and_then(|arc| decode(arc.start))
        .ok_or("empty circle group")?;
    for seed in seeds {
        if visited.contains(&seed) {
            continue;
        }
        let start = arcs[seed].start;
        let mut current = seed;
        let mut total = 0.;
        let end = loop {
            if !visited.insert(current) {
                return Err("repeated circular boundary".into());
            }
            let arc = arcs[current];
            total += arc.sweep;
            if total > TAU + tol {
                return Err("circular boundary winds more than once".into());
            }
            if arc.end == start {
                break arc.end;
            }
            let Some(&next) = outgoing.get(&arc.end) else {
                break arc.end;
            };
            current = next;
        };
        let endpoints_mm = if start == end {
            if (total - TAU).abs() > tol {
                return Err("incomplete circular cycle".into());
            }
            None
        } else {
            let expected = positive_angle(
                radial(decode(start).ok_or("invalid arc start")?),
                radial(decode(end).ok_or("invalid arc end")?),
                normal,
            );
            if (total - expected).abs() > tol {
                return Err("inconsistent directed circular arc coverage".into());
            }
            Some([start, end])
        };
        if endpoints_mm.is_none() {
            // A certified complete revolution has one interval, even if
            // its accumulated span sweep exceeds TAU by roundoff.
            intervals.push((0.0, TAU));
        } else {
            let mut from = positive_angle(
                radial(angle_origin),
                radial(decode(start).ok_or("invalid arc start")?),
                normal,
            );
            if TAU - from <= tol {
                from = 0.0;
            }
            let to = from + total;
            intervals.push((from, to.min(TAU)));
            if to > TAU + tol {
                intervals.push((0.0, to - TAU));
            }
        }
        result.push(endpoints_mm);
    }
    intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
    if intervals.windows(2).any(|pair| pair[0].1 > pair[1].0 + tol) {
        return Err("overlapping circular boundary coverage".into());
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::super::point_key;
    use super::*;
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI};
    const Z: [i64; 3] = [0, 0, 10_000_000_000];
    const R: i64 = 10_000_000;

    fn plane() -> OrientedSupportKey {
        OrientedSupportKey::Plane {
            normal: Z,
            offset_mm: 0,
            outward_positive: true,
        }
    }
    fn sphere() -> OrientedSupportKey {
        OrientedSupportKey::Sphere {
            center_mm: [0; 3],
            radius_mm: R,
            radial_outward: true,
        }
    }
    fn unit_circle() -> brep::CircleSupport {
        brep::CircleSupport {
            center_mm: [0.; 3],
            normal: [0., 0., 1.],
            x_direction: [1., 0., 0.],
            radius_mm: 1.,
        }
    }
    fn edge(
        start: PointKey,
        end: PointKey,
        parameter_forward: bool,
        curve: CurveSupport,
    ) -> EdgeGeometry {
        EdgeGeometry {
            start,
            end,
            start_vertex: 1,
            end_vertex: 2,
            parameter_forward,
            curve_id: 10,
            support: curve,
        }
    }
    fn arc(start: PointKey, end: PointKey, forward: bool) -> anyhow::Result<(CircleLocusKey, Arc)> {
        let curve = CurveSupport::Circle(unit_circle());
        let edge = edge(start, end, forward, curve.clone());
        let loci = [support_locus_key(&plane()), support_locus_key(&sphere())];
        let Some(BoundedEdge::Circle(locus, arc)) =
            classify(&edge, &curve, &loci, &mut CellWitnesses::default())
        else {
            anyhow::bail!("expected certified circular arc")
        };
        Ok((locus, arc))
    }

    #[test]
    fn complementary_arcs_are_distinct_and_recharting_preserves_arc() -> anyhow::Result<()> {
        let supports = [plane(), sphere()];
        let start = [R, 0, 0];
        let end = [0, R, 0];
        let (locus, short) = arc(start, end, true)?;
        let (_, long) = arc(start, end, false)?;
        let (_, reversed) = arc(end, start, false)?;
        let signature = |arc| normalize_arcs(&supports, &locus, &[arc]);
        assert_eq!(signature(short), signature(reversed));
        assert_ne!(signature(short), signature(long));
        assert!((short.sweep - FRAC_PI_2).abs() < TOL);
        assert!((long.sweep - (TAU - FRAC_PI_2)).abs() < TOL);
        Ok(())
    }

    #[test]
    fn directed_circle_chain_ignores_partitions_but_rejects_overlaps() -> anyhow::Result<()> {
        let supports = [plane(), sphere()];
        let start = [R, 0, 0];
        let middle = [0, R, 0];
        let end = [-R, 0, 0];
        let (locus, half) = arc(start, end, true)?;
        let (_, first) = arc(start, middle, true)?;
        let (_, second) = arc(middle, end, true)?;
        let (_, other_half) = arc(end, start, true)?;
        assert_eq!(
            normalize_arcs(&supports, &locus, &[half]),
            normalize_arcs(&supports, &locus, &[second, first])
        );
        let cycle =
            normalize_arcs(&supports, &locus, &[half, other_half]).map_err(anyhow::Error::msg)?;
        let rounded = Arc {
            sweep: 4.0_f64.mul_add(f64::EPSILON, other_half.sweep),
            ..other_half
        };
        assert_eq!(
            normalize_arcs(&supports, &locus, &[half, rounded]).map_err(anyhow::Error::msg)?,
            cycle
        );

        assert!(matches!(
            cycle.as_slice(),
            [BoundaryChainKey::Circular {
                endpoints_mm: None,
                ..
            }]
        ));
        assert!(normalize_arcs(&supports, &locus, &[half, first, second]).is_err());
        let (_, complement) = arc(middle, start, true)?;
        assert!(normalize_arcs(&supports, &locus, &[half, complement]).is_err());
        Ok(())
    }

    #[test]
    fn sphere_and_torus_circle_incidence_requires_whole_carrier() {
        let sphere = support_locus_key(&sphere());
        assert!(circle_on_curved_support(unit_circle(), &sphere));
        assert!(!circle_on_curved_support(
            brep::CircleSupport {
                center_mm: [0., 0., 0.01],
                ..unit_circle()
            },
            &sphere
        ));
        let torus = SupportLocusKey::Torus {
            center_mm: [0; 3],
            axis: Z,
            major_radius_mm: 3 * R,
            minor_radius_mm: R,
        };
        let parallel = brep::CircleSupport {
            radius_mm: 3.,
            center_mm: [0., 0., 1.],
            ..unit_circle()
        };
        let meridian = brep::CircleSupport {
            center_mm: [3., 0., 0.],
            normal: [0., 1., 0.],
            ..unit_circle()
        };
        assert!(circle_on_curved_support(parallel, &torus));
        assert!(circle_on_curved_support(meridian, &torus));
        assert!(!circle_on_curved_support(
            brep::CircleSupport {
                radius_mm: 0.99,
                ..meridian
            },
            &torus
        ));
        // Two distinct circular components of a torus-plane intersection:
        // identical oriented carrier pair cannot identify their loci.
        let inner = brep::CircleSupport {
            radius_mm: 2.,
            ..unit_circle()
        };
        let outer = brep::CircleSupport {
            radius_mm: 4.,
            ..unit_circle()
        };
        assert!(circle_on_curved_support(inner, &torus));
        assert!(circle_on_curved_support(outer, &torus));
        let loci = [support_locus_key(&plane()), torus];
        let classify_circle = |c: brep::CircleSupport| {
            let curve = CurveSupport::Circle(c);
            let mut edge = edge(
                point_key([c.radius_mm, 0., 0.]).unwrap_or_default(),
                point_key([c.radius_mm, 0., 0.]).unwrap_or_default(),
                true,
                curve.clone(),
            );
            edge.end_vertex = edge.start_vertex;
            match classify(&edge, &curve, &loci, &mut CellWitnesses::default()) {
                Some(BoundedEdge::Circle(key, _)) => Some(key),
                _ => None,
            }
        };
        let first = classify_circle(inner);
        let second = classify_circle(outer);
        assert!(first.is_some() && second.is_some());
        assert_ne!(first, second);
    }

    fn oblique_supports() -> [OrientedSupportKey; 2] {
        [
            OrientedSupportKey::Plane {
                normal: [-7_071_067_812, 0, 7_071_067_812],
                offset_mm: 0,
                outward_positive: true,
            },
            OrientedSupportKey::Cylinder {
                axis: Z,
                closest_mm: [0; 3],
                radius_mm: R,
                radial_outward: true,
            },
        ]
    }

    fn ellipse_arc(lo: f64, hi: f64) -> BSplineSupport {
        let point = |a: f64, r: f64| {
            let x = r * a.cos();
            [x, r * a.sin(), x]
        };
        let w = ((hi - lo) * 0.5).cos();
        BSplineSupport {
            degree: 2,
            control_points_mm: vec![
                point(lo, 1.),
                point(lo.midpoint(hi), w.recip()),
                point(hi, 1.),
            ],
            knots: vec![0., 0., 0., 1., 1., 1.],
            weights: Some(vec![1., w, 1.]),
        }
    }

    fn oblique_arc(lo: f64, hi: f64) -> anyhow::Result<Arc> {
        let spline = ellipse_arc(lo, hi);
        let start = point_key(
            *spline
                .control_points_mm
                .first()
                .ok_or_else(|| anyhow::anyhow!("no start"))?,
        )
        .ok_or_else(|| anyhow::anyhow!("invalid start"))?;
        let end = point_key(
            *spline
                .control_points_mm
                .last()
                .ok_or_else(|| anyhow::anyhow!("no end"))?,
        )
        .ok_or_else(|| anyhow::anyhow!("invalid end"))?;
        let curve = CurveSupport::BSpline(spline);
        let e = edge(start, end, true, curve.clone());
        let Some(BoundedEdge::CylinderSection(arc)) = classify(
            &e,
            &curve,
            &oblique_supports().each_ref().map(support_locus_key),
            &mut CellWitnesses::default(),
        ) else {
            anyhow::bail!("exact oblique ellipse rejected");
        };
        Ok(arc)
    }

    #[test]
    fn oblique_ellipse_coverage_is_partition_independent() -> anyhow::Result<()> {
        let supports = oblique_supports();
        let whole = oblique_arc(0., FRAC_PI_2)?;
        let first = oblique_arc(0., FRAC_PI_2 * 0.5)?;
        let second = oblique_arc(FRAC_PI_2 * 0.5, FRAC_PI_2)?;
        let reversed = oblique_arc(FRAC_PI_2, 0.)?;
        let signature = normalize_cylinder_sections(&supports, &[whole]);
        assert_eq!(
            signature,
            normalize_cylinder_sections(&supports, &[second, first])
        );
        assert_eq!(
            signature,
            normalize_cylinder_sections(&supports, &[reversed])
        );
        assert!(signature.is_ok());
        assert!(normalize_cylinder_sections(&supports, &[whole, first, second]).is_err());
        assert!(matches!(
            signature.map_err(anyhow::Error::msg)?.as_slice(),
            [BoundaryChainKey::CylinderSection {
                endpoints_mm: Some([[R, 0, R], [0, R, 0]]),
                ..
            }]
        ));
        let rest = [
            oblique_arc(FRAC_PI_2, PI)?,
            oblique_arc(PI, PI + FRAC_PI_2)?,
            oblique_arc(PI + FRAC_PI_2, TAU)?,
        ];
        assert!(matches!(
            normalize_cylinder_sections(&supports, &[whole, rest[0], rest[1], rest[2]])
                .map_err(anyhow::Error::msg)?
                .as_slice(),
            [BoundaryChainKey::CylinderSection {
                endpoints_mm: None,
                ..
            }]
        ));
        Ok(())
    }

    #[test]
    fn declared_self_intersection_does_not_override_complete_ellipse_proof() -> anyhow::Result<()> {
        let input = "ISO-10303-21;HEADER;FILE_DESCRIPTION(('x'),'1');FILE_NAME('a','b',(''),(''),'x','y','');FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;
        #1=CARTESIAN_POINT('',(1.,0.,1.));
        #2=CARTESIAN_POINT('',(1.,1.,1.));
        #3=CARTESIAN_POINT('',(0.,1.,0.));
        #4=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#1,#2,#3),.UNSPECIFIED.,.F.,.T.)B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,0.7071067811865476,1.))REPRESENTATION_ITEM(''));
        #5=INTERSECTION_CURVE('',#4,(),.CURVE_3D.);
        ENDSEC;END-ISO-10303-21;";
        let exchange = step_p21::parser::parse(input)?;
        let entities = &exchange.data[0].entities;
        let index = crate::step_graph::build_index(entities);
        let ordinary = brep::curve_support(5, entities, &index);
        assert!(matches!(ordinary, CurveSupport::Other { .. }));
        let explicit = brep::curve_geometry_for_incidence_proof(5, entities, &index);
        let e = edge([R, 0, R], [0, R, 0], true, ordinary);
        let loci = oblique_supports().each_ref().map(support_locus_key);
        assert!(matches!(
            certify_flagged_cylinder_section(&e, &explicit, &loci),
            Some(BoundedEdge::CylinderSection(_))
        ));
        let CurveSupport::BSpline(mut wrong) = explicit else {
            anyhow::bail!("expected rational input");
        };
        wrong.control_points_mm[1][2] += 0.01;
        assert!(
            certify_flagged_cylinder_section(&e, &CurveSupport::BSpline(wrong), &loci).is_none()
        );
        Ok(())
    }

    #[test]
    fn oblique_ellipse_requires_entire_plane_and_cylinder_incidence() {
        let supports = oblique_supports();
        let loci = supports.each_ref().map(support_locus_key);
        let spline = ellipse_arc(0., FRAC_PI_2);
        let e = edge(
            [R, 0, R],
            [0, R, 0],
            true,
            CurveSupport::BSpline(spline.clone()),
        );
        let mut off_plane = spline.clone();
        off_plane.control_points_mm[1][2] += 0.001;
        assert!(
            classify(
                &e,
                &CurveSupport::BSpline(off_plane),
                &loci,
                &mut CellWitnesses::default()
            )
            .is_none()
        );
        let mut off_cylinder = spline.clone();
        off_cylinder.control_points_mm[1][0] += 0.001;
        off_cylinder.control_points_mm[1][2] += 0.001;
        assert!(
            classify(
                &e,
                &CurveSupport::BSpline(off_cylinder),
                &loci,
                &mut CellWitnesses::default()
            )
            .is_none()
        );
        let mut wrong = e;
        wrong.end = [-R, 0, -R];
        assert!(
            classify(
                &wrong,
                &CurveSupport::BSpline(spline),
                &loci,
                &mut CellWitnesses::default()
            )
            .is_none()
        );
    }

    #[test]
    fn rational_arc_certificate_rejects_off_carrier_and_wrong_endpoints() -> anyhow::Result<()> {
        let supports = [plane(), sphere()];
        let loci = supports.each_ref().map(support_locus_key);
        let curve = BSplineSupport {
            degree: 2,
            control_points_mm: vec![[1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
            knots: vec![0., 0., 0., 1., 1., 1.],
            weights: Some(vec![1., FRAC_1_SQRT_2, 1.]),
        };
        let exact = CurveSupport::BSpline(curve.clone());
        let edge = edge([R, 0, 0], [0, R, 0], true, exact.clone());
        let Some(BoundedEdge::Circle(locus, certified)) =
            classify(&edge, &exact, &loci, &mut CellWitnesses::default())
        else {
            anyhow::bail!("exact rational circle refused")
        };
        let (reference_locus, reference) = arc([R, 0, 0], [0, R, 0], true)?;
        assert_eq!(locus, reference_locus);
        assert_eq!(
            normalize_arcs(&supports, &locus, &[certified]),
            normalize_arcs(&supports, &locus, &[reference])
        );
        let mut bumped = curve;
        bumped.control_points_mm[1][0] += 0.001;
        assert!(
            classify(
                &edge,
                &CurveSupport::BSpline(bumped),
                &loci,
                &mut CellWitnesses::default()
            )
            .is_none()
        );
        let mut wrong = edge;
        wrong.end = [-R, 0, 0];
        assert!(classify(&wrong, &exact, &loci, &mut CellWitnesses::default()).is_none());
        assert!((certified.sweep - PI / 2.).abs() < TOL);
        Ok(())
    }
}
