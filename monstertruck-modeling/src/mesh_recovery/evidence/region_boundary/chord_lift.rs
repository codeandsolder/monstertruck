//! Exact single-carrier paths for boundary-constrained approximations.
//!
//! A lifted chord is not an inferred original intersection. It replaces a
//! source segment by a rational path on one recognized carrier while the
//! adjoining unresolved facet may use a separately bounded approximation.

use super::{SurfaceCarrier, SurfaceRegionBudget};
use crate::{BsplineCurve, InnerSpace, KnotVector, NurbsCurve, Point3, Vector3, Vector4};
use monstertruck_geometry::analytic_line_segment::bound_nurbs_segment_locus;

/// A source segment lifted to an independently fitted analytic carrier.
#[derive(Clone, Debug)]
pub struct LiftedSourceChord {
    /// Rational curve on the selected carrier, in source endpoint order.
    pub curve: NurbsCurve<Vector4>,
    /// Exact projected junction positions used by the curve.
    pub endpoints: [Point3; 2],
    /// Maximum displacement of either original endpoint.
    pub vertex_displacement: f64,
    /// Continuous two-sided deviation from the original finite segment.
    pub chord_deviation_upper_bound: f64,
}

fn finite(p: Point3) -> bool { [p.x, p.y, p.z].into_iter().all(f64::is_finite) }

fn unit(axis: Vector3) -> bool { (axis.magnitude2() - 1.0).abs() <= 64.0 * f64::EPSILON }

fn homogeneous(point: Point3, w: f64) -> Vector4 {
    Vector4::new(point.x * w, point.y * w, point.z * w, w)
}

fn cartesian(point: Vector4) -> Point3 {
    Point3::new(point.x / point.w, point.y / point.w, point.z / point.w)
}

fn bezier(controls: Vec<Vector4>) -> Option<NurbsCurve<Vector4>> {
    if controls.len() < 2
        || controls
            .iter()
            .any(|p| p.w <= 0.0 || ![p.x, p.y, p.z, p.w].into_iter().all(f64::is_finite))
    {
        return None;
    }
    Some(NurbsCurve::new(
        BsplineCurve::new(KnotVector::bezier_knot(controls.len() - 1), controls).ok()?,
    ))
}

#[derive(Clone, Copy)]
struct WeightedDirection {
    numerator: Vector3,
    weight: f64,
}

// Short unit-circle arc in an arbitrary plane through the origin. The
// homogeneous middle numerator is the unit angular bisector, not the
// Euclidean tangent intersection. Both coincident directions and oblique
// great circles are handled by the same formula.
fn direction_arc(a: Vector3, b: Vector3) -> Option<[WeightedDirection; 3]> {
    let bisector = a + b;
    let weight = bisector.magnitude() * 0.5;
    if !weight.is_finite() || weight <= 128.0 * f64::EPSILON {
        return None; // Antipodal: samples do not determine a short branch.
    }
    Some([
        WeightedDirection {
            numerator: a,
            weight: 1.0,
        },
        WeightedDirection {
            numerator: bisector / (2.0 * weight),
            weight,
        },
        WeightedDirection {
            numerator: b,
            weight: 1.0,
        },
    ])
}

#[derive(Clone, Copy)]
struct AxialPoint {
    axial: f64,
    radius: f64,
    direction: Option<Vector3>,
}

fn axial_point(origin: Point3, axis: Vector3, point: Point3) -> Option<AxialPoint> {
    let delta = point - origin;
    let axial = delta.dot(axis);
    let radial = delta - axis * axial;
    let radius = radial.magnitude();
    if !axial.is_finite() || !radius.is_finite() {
        return None;
    }
    let direction = (radius > 0.0).then(|| radial / radius);
    Some(AxialPoint {
        axial,
        radius,
        direction,
    })
}

fn azimuth(a: AxialPoint, b: AxialPoint) -> Option<[WeightedDirection; 3]> {
    // A point at an actual rotational pole can use the other endpoint's
    // azimuth. An off-axis positive radius must have its own direction.
    let first = a
        .direction
        .or_else(|| (a.radius == 0.0).then_some(b.direction).flatten())?;
    let last = b
        .direction
        .or_else(|| (b.radius == 0.0).then_some(a.direction).flatten())?;
    direction_arc(first, last)
}

fn plane(origin: Point3, normal: Vector3, source: [Point3; 2]) -> Option<Vec<Vector4>> {
    if !finite(origin) || !unit(normal) {
        return None;
    }
    Some(
        source
            .into_iter()
            .map(|p| homogeneous(p - normal * (p - origin).dot(normal), 1.0))
            .collect(),
    )
}

fn sphere(center: Point3, radius: f64, source: [Point3; 2]) -> Option<Vec<Vector4>> {
    if !finite(center) || !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let [a, b] = source.map(|p| p - center);
    if a.magnitude() <= 0.0 || b.magnitude() <= 0.0 {
        return None;
    }
    let arc = direction_arc(a.normalize(), b.normalize())?;
    Some(
        arc.into_iter()
            .map(|p| {
                let n = p.numerator * radius;
                homogeneous(center, p.weight) + Vector4::new(n.x, n.y, n.z, 0.0)
            })
            .collect(),
    )
}

fn project_linear(p: AxialPoint, reference_radius: f64, slope: f64) -> Option<AxialPoint> {
    let norm = 1.0_f64.hypot(slope);
    let residual = (p.radius - slope.mul_add(p.axial, reference_radius)) / norm;
    let axial = (slope / norm).mul_add(residual, p.axial);
    let radius = slope.mul_add(axial, reference_radius);
    if !axial.is_finite() || !radius.is_finite() {
        return None;
    }
    if radius < 0.0 {
        // Closest point of the positive-radius half-line is the apex.
        if slope == 0.0 {
            return None;
        }
        let axial = -reference_radius / slope;
        return axial.is_finite().then_some(AxialPoint {
            axial,
            radius: 0.0,
            direction: p.direction,
        });
    }
    Some(AxialPoint {
        axial,
        radius,
        direction: p.direction,
    })
}

fn linear(
    origin: Point3,
    axis: Vector3,
    reference_radius: f64,
    slope: f64,
    source: [Point3; 2],
) -> Option<Vec<Vector4>> {
    if !finite(origin) || !unit(axis) || !reference_radius.is_finite() || !slope.is_finite() {
        return None;
    }
    let a = project_linear(
        axial_point(origin, axis, source[0])?,
        reference_radius,
        slope,
    )?;
    let b = project_linear(
        axial_point(origin, axis, source[1])?,
        reference_radius,
        slope,
    )?;
    let arc = azimuth(a, b)?;
    // Exact product of the quadratic rational azimuth with a linear
    // meridian. Use Bernstein product factors C(2,i)C(1,j)/C(3,i+j).
    let mut controls = vec![Vector4::new(0.0, 0.0, 0.0, 0.0); 4];
    for (i, p) in arc.into_iter().enumerate() {
        for (j, m) in [a, b].into_iter().enumerate() {
            let factor = [1.0, 2.0, 1.0][i] / [1.0, 3.0, 3.0, 1.0][i + j];
            let offset = axis * (m.axial * p.weight) + p.numerator * m.radius;
            controls[i + j] += (homogeneous(origin, p.weight)
                + Vector4::new(offset.x, offset.y, offset.z, 0.0))
                * factor;
        }
    }
    Some(controls)
}

fn meridian(point: AxialPoint, center_radius: f64) -> Option<Vector3> {
    let p = Vector3::new(point.radius - center_radius, point.axial, 0.0);
    let norm = p.magnitude();
    (norm.is_finite() && norm > 0.0).then(|| p / norm)
}

// Minimum cosine on the selected short angular interval, including
// stationary points. This matters for spindle tori: an admissible pair
// of endpoints alone does not prove the whole arc stays on the positive
// radial sheet.
fn min_cosine(a: Vector3, b: Vector3) -> f64 {
    use std::f64::consts::PI;
    let first = a.y.atan2(a.x);
    let sweep = a.x.mul_add(b.y, -a.y * b.x).atan2(a.dot(b));
    let last = first + sweep;
    let lo = first.min(last);
    let hi = first.max(last);
    if (lo..=hi).contains(&PI) || (lo..=hi).contains(&-PI) {
        -1.0
    } else {
        a.x.min(b.x)
    }
}

fn torus(
    center: Point3,
    axis: Vector3,
    center_radius: f64,
    profile_radius: f64,
    source: [Point3; 2],
) -> Option<Vec<Vector4>> {
    if !finite(center)
        || !unit(axis)
        || !center_radius.is_finite()
        || center_radius <= 0.0
        || !profile_radius.is_finite()
        || profile_radius <= 0.0
    {
        return None;
    }
    let a = axial_point(center, axis, source[0])?;
    let b = axial_point(center, axis, source[1])?;
    let u = meridian(a, center_radius)?;
    let v = meridian(b, center_radius)?;
    let profile = direction_arc(u, v)?;
    let radial_min = profile_radius.mul_add(min_cosine(u, v), center_radius);
    if radial_min <= 128.0 * f64::EPSILON * center_radius.max(profile_radius) {
        return None; // No regular positive-radial chart over this branch.
    }
    let azimuth = direction_arc(a.direction?, b.direction?)?;
    // Tensor product evaluated along a common rational parameter. Both
    // angles follow short arcs exactly, so the resulting quartic lies on
    // the toroidal carrier even when neither angle is constant.
    let mut controls = vec![Vector4::new(0.0, 0.0, 0.0, 0.0); 5];
    for (i, p) in azimuth.into_iter().enumerate() {
        for (j, m) in profile.into_iter().enumerate() {
            let factor = [1.0, 2.0, 1.0][i] * [1.0, 2.0, 1.0][j] / [1.0, 4.0, 6.0, 4.0, 1.0][i + j];
            let weight = p.weight * m.weight;
            let radial = profile_radius.mul_add(m.numerator.x, center_radius * m.weight);
            let offset = p.numerator * radial + axis * (profile_radius * m.numerator.y * p.weight);
            controls[i + j] += (homogeneous(center, weight)
                + Vector4::new(offset.x, offset.y, offset.z, 0.0))
                * factor;
        }
    }
    Some(controls)
}

impl SurfaceCarrier {
    /// Lift one finite source segment onto this carrier using an exact
    /// rational chart path. The source endpoints are projected first;
    /// each must respect the vertex budget independently of the continuous
    /// chord-deviation budget. Antipodal and singular chart branches that
    /// the two samples cannot resolve are rejected.
    ///
    /// This is explicit approximation evidence, not a claim that the
    /// original CAD edge was the returned path.
    #[must_use]
    pub fn lift_source_chord(
        self,
        source: [Point3; 2],
        budget: SurfaceRegionBudget,
    ) -> Option<LiftedSourceChord> {
        if !source.into_iter().all(finite)
            || !budget.vertex.is_finite()
            || budget.vertex < 0.0
            || !budget.chord.is_finite()
            || budget.chord <= 0.0
            || (source[1] - source[0]).magnitude() <= 0.0
        {
            return None;
        }
        let controls = match self {
            Self::Plane { origin, normal } => plane(origin, normal, source)?,
            Self::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } => linear(origin, axis, radius_at_origin, slope, source)?,
            Self::CircularRevolution {
                axial_center,
                radial_center: 0.0,
                profile_radius,
                ..
            } => sphere(axial_center, profile_radius, source)?,
            Self::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                profile_radius,
            } => torus(axial_center, axis, radial_center, profile_radius, source)?,
        };
        let endpoints = [cartesian(*controls.first()?), cartesian(*controls.last()?)];
        let vertex_displacement = (endpoints[0] - source[0])
            .magnitude()
            .max((endpoints[1] - source[1]).magnitude());
        if !vertex_displacement.is_finite()
            || vertex_displacement > budget.vertex
            || (endpoints[1] - endpoints[0]).magnitude() <= 0.0
        {
            return None;
        }
        let curve = bezier(controls)?;
        let chord_deviation_upper_bound =
            bound_nurbs_segment_locus(&curve, source[0], source[1], budget.chord)?;
        (chord_deviation_upper_bound <= budget.chord).then_some(LiftedSourceChord {
            curve,
            endpoints,
            vertex_displacement,
            chord_deviation_upper_bound,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ParametricCurve;

    fn check(carrier: SurfaceCarrier, source: [Point3; 2], max_degree: usize) {
        let result = carrier
            .lift_source_chord(
                source,
                SurfaceRegionBudget {
                    vertex: 0.001,
                    chord: 0.5,
                },
            )
            .expect("regular lifted chord");
        assert!(result.curve.non_rationalized().degree() <= max_degree);
        assert!(result.vertex_displacement <= 0.001);
        for i in 0..=100 {
            let p = result.curve.evaluate(f64::from(i) / 100.0);
            let (_, residual) = carrier.signed_constraint(p).expect("regular carrier");
            assert!(residual.abs() < 1.0e-11, "{residual}");
            let direction = source[1] - source[0];
            let t = ((p - source[0]).dot(direction) / direction.magnitude2()).clamp(0.0, 1.0);
            assert!(
                (p - (source[0] + direction * t)).magnitude()
                    <= result.chord_deviation_upper_bound + 1.0e-12
            );
        }
    }

    #[test]
    fn single_carrier_paths_are_exact_on_all_regular_families() {
        let origin = Point3::new(3.0, -5.0, 7.0);
        let axis = Vector3::unit_z();
        let x = Vector3::unit_x();
        let y = Vector3::unit_y();
        check(
            SurfaceCarrier::Plane {
                origin,
                normal: axis,
            },
            [origin + x, origin + y + axis * 0.0001],
            1,
        );
        check(
            SurfaceCarrier::CircularRevolution {
                axial_center: origin,
                axis,
                radial_center: 0.0,
                profile_radius: 2.0,
            },
            [
                origin + x * 2.0,
                origin + (x + y + axis).normalize() * 2.0001,
            ],
            2,
        );
        for slope in [0.0, 0.4, -0.3] {
            check(
                SurfaceCarrier::LinearRevolution {
                    origin,
                    axis,
                    radius_at_origin: 2.0,
                    slope,
                },
                [
                    origin + x * 2.0001,
                    origin + axis * 0.5 + (x + y).normalize() * slope.mul_add(0.5, 2.0),
                ],
                3,
            );
        }
        for center_radius in [3.0, 0.8] {
            let point = |u: f64, v: f64| {
                let r = v.cos() + center_radius;
                origin + (x * u.cos() + y * u.sin()) * r + axis * v.sin()
            };
            check(
                SurfaceCarrier::CircularRevolution {
                    axial_center: origin,
                    axis,
                    radial_center: center_radius,
                    profile_radius: 1.0,
                },
                [point(0.2, -0.3), point(0.5, 0.4)],
                4,
            );
        }
    }

    #[test]
    fn independent_budgets_and_ambiguous_branches_are_rejected() {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: origin,
            axis: Vector3::unit_z(),
            radial_center: 0.0,
            profile_radius: 1.0,
        };
        let a = Point3::new(1.0, 0.0, 0.0);
        let b = Point3::new(0.0, 1.0, 0.0);
        assert!(
            carrier
                .lift_source_chord(
                    [a, b],
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        chord: 0.2
                    }
                )
                .is_none()
        );
        assert!(
            carrier
                .lift_source_chord(
                    [a, b],
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        chord: 0.3
                    }
                )
                .is_some()
        );
        assert!(
            carrier
                .lift_source_chord(
                    [a, Point3::new(-1.0, 0.0, 0.0)],
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        chord: 2.0
                    }
                )
                .is_none()
        );
        assert!(
            carrier
                .lift_source_chord(
                    [a, Point3::new(0.0, 1.1, 0.0)],
                    SurfaceRegionBudget {
                        vertex: 0.01,
                        chord: 1.0
                    }
                )
                .is_none()
        );
        // Both endpoints lie on the positive-radius sheet of a spindle,
        // but the shorter meridian arc through cos(v)=-1 crosses its axis.
        let spindle = SurfaceCarrier::CircularRevolution {
            axial_center: origin,
            axis: Vector3::unit_z(),
            radial_center: 0.5,
            profile_radius: 1.0,
        };
        let angle = 1.9_f64;
        let endpoints = [
            Point3::new(0.5 + angle.cos(), 0.0, angle.sin()),
            Point3::new(0.5 + angle.cos(), 0.0, -angle.sin()),
        ];
        assert!(
            spindle
                .lift_source_chord(
                    endpoints,
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        chord: 3.0
                    }
                )
                .is_none()
        );
    }
}
