//! Boundary-constrained approximations of unresolved source facets.
//!
//! This is not spline-carrier recognition: sparse mesh samples generally
//! do not determine the original freeform surface. Given certified shared
//! boundary edges, construct an explicitly bounded replacement for a
//! convex source facet. The shared edge itself is never approximated again.

mod apex;

use crate::{
    BsplineCurve, BsplineSurface, Edge, EuclideanSpace, Face, InnerSpace, KnotVector, NurbsSurface,
    Point3, Surface, Vector3, Vector4, Wire, builder,
};
use monstertruck_geometry::{
    prelude::TryIntoHomogeneousBsplineCurve, proof_bezier_spans, subdivide_homogeneous_bezier,
};

/// A piecewise rational approximation of a convex source facet.
///
/// The source facet's boundary is retained using the supplied native edge
/// identities. Internal radial edges are shared between the fan faces.
#[derive(Debug)]
pub struct BoundaryFacetApproximation {
    /// Oriented native fan faces, in source boundary order.
    pub faces: Vec<Face>,
    /// Continuous two-sided distance bound to the source polygon's
    /// reference centroid fan, in source coordinate units.
    pub deviation_upper_bound: f64,
    /// Maximum bounded discrepancy between a radial seam and the exact
    /// curve endpoint fan, including the numerical certificate margin.
    pub seam_deviation_upper_bound: f64,
}

/// Why a boundary-constrained facet cannot be certified.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum BoundaryFacetError {
    /// Invalid input, nonconvex/degenerate source or displaced boundary.
    #[error("invalid source facet or native boundary")]
    InvalidBoundary,
    /// The exact native curve lacks a positive-weight homogeneous net.
    #[error("unsupported exact boundary curve")]
    UnsupportedCurve,
    /// The curve leaves its angular wedge or exceeds the approximation budget.
    #[error("boundary approximation or orientation could not be certified")]
    Uncertified,
    /// The native kernel rejected the constructed fan face.
    #[error("invalid native fan face")]
    NativeFace,
}

fn finite(point: Point3) -> bool { [point.x, point.y, point.z].into_iter().all(f64::is_finite) }

fn centroid(points: &[Point3]) -> Point3 {
    let first = points[0];
    first
        + points
            .iter()
            .map(|&p| p - first)
            .fold(Vector3::new(0.0, 0.0, 0.0), |sum, p| sum + p)
            / crate::usize_to_f64(points.len())
}

fn directed_pairs<T>(values: &[T]) -> impl Iterator<Item = (&T, &T)> {
    values.iter().zip(values.iter().cycle().skip(1))
}

fn cross(normal: Vector3, a: Vector3, b: Vector3) -> f64 { normal.dot(a.cross(b)) }

// Require a convex polygon and an interior apex in its projected plane.
// With strict adjacent ray order this partitions the plane into disjoint
// angular wedges, regardless of the actual edge parameter speed.
fn convex_fan(points: &[Point3], apex: Point3, normal: Vector3) -> bool {
    directed_pairs(points).all(|(&a, &b)| {
        let edge = b - a;
        let scale = edge.magnitude();
        let allowance = 64.0 * f64::EPSILON * scale;
        scale > 0.0
            && cross(normal, edge, apex - a) > allowance * (apex - a).magnitude()
            && points
                .iter()
                .all(|&p| cross(normal, edge, p - a) >= -allowance * (p - a).magnitude())
    })
}

fn chord_distance(point: Point3, a: Point3, b: Point3) -> f64 {
    let delta = b - a;
    let t = ((point - a).dot(delta) / delta.magnitude2()).clamp(0.0, 1.0);
    (point - (a + delta * t)).magnitude()
}

fn rational_point(control: Vector4) -> Point3 {
    Point3::new(
        control.x / control.w,
        control.y / control.w,
        control.z / control.w,
    )
}

fn positive_weights(controls: &[Vector4]) -> bool {
    controls
        .iter()
        .all(|&p| p.w > 0.0 && [p.x, p.y, p.z, p.w].into_iter().all(f64::is_finite))
}

struct Wedge {
    apex: Point3,
    normal: Vector3,
    source: [Point3; 2],
    rays: [Vector3; 2],
    budget: f64,
}

// The native topology uses shared vertex rays. A fitted analytic curve
// can have a distinct endpoint within the kernel's incidence tolerance
// (e.g. independently fitted near-tangent carriers). Prove injectivity in
// the actual curve's angular wedge AND confinement to the shared wedge
// up to an explicitly measured endpoint discrepancy.
struct AngularWedge {
    half_plane: Vector3,
    shared_rays: [Vector3; 2],
    seam_gaps: [f64; 2],
}

struct FacetCurveCertificate {
    deviation: f64,
    seam_deviation: f64,
}

impl Wedge {
    // All projected control points must remain in the wedge bounded by
    // the two endpoint rays. The rational angular derivative has the
    // sign of n·(P×P'), where P is the homogeneous numerator measured
    // relative to the apex. The derivative of the denominator cancels.
    // A positive bound on every control/derivative pair certifies the
    // entire span by the Bernstein convex-hull property.
    fn monotone(&self, controls: &[Vector4], allowance: f64, angular: &AngularWedge) -> bool {
        let max_weight = controls.iter().map(|p| p.w).fold(0.0_f64, f64::max);
        let local = controls
            .iter()
            .map(|p| {
                Vector3::new(
                    (-self.apex.x).mul_add(p.w, p.x),
                    (-self.apex.y).mul_add(p.w, p.y),
                    (-self.apex.z).mul_add(p.w, p.z),
                ) / max_weight
            })
            .collect::<Vec<_>>();
        let in_wedge = controls.iter().all(|&control| {
            // Euclidean control hulls certify actual distances. In
            // particular, never apply an unweighted homogeneous error
            // margin as if it were a spatial seam-distance bound.
            let p = rational_point(control) - self.apex;
            let roundoff = 2.0 * allowance;
            let [a, b] = angular.shared_rays;
            let [gap_a, gap_b] = angular.seam_gaps;
            p.dot(angular.half_plane) > roundoff
                && cross(self.normal, a, p) >= -(gap_a + roundoff) * a.magnitude()
                && cross(self.normal, p, b) >= -(gap_b + roundoff) * b.magnitude()
        });
        in_wedge
            && local.windows(2).all(|pair| {
                let derivative = pair[1] - pair[0];
                let derivative_norm = derivative.magnitude();
                local.iter().all(|&p| {
                    let norm = p.magnitude();
                    let uncertainty =
                        2.0_f64.mul_add(allowance, 2.0_f64.mul_add(norm, derivative_norm));
                    let error = (64.0 * f64::EPSILON * norm)
                        .mul_add(derivative_norm, allowance * uncertainty);
                    cross(self.normal, p, derivative) > error
                })
            })
    }

    fn certify(
        &self,
        controls: &[[f64; 4]],
        allowance: f64,
        depth: usize,
        remaining: &mut usize,
        angular: &AngularWedge,
    ) -> Option<f64> {
        *remaining = remaining.checked_sub(1)?;
        let controls_vec = controls
            .iter()
            .map(|&[x, y, z, w]| Vector4::new(x, y, z, w))
            .collect::<Vec<_>>();
        if controls_vec.len() < 2 || !positive_weights(&controls_vec) {
            return None;
        }
        // Positive rational weights keep the entire curve in the convex
        // hull of Euclidean controls. Distance to a closed line segment
        // is convex, so the maximum control distance is a valid bound.
        let distances = controls_vec
            .iter()
            .map(|&p| chord_distance(rational_point(p), self.source[0], self.source[1]))
            .collect::<Vec<_>>();
        let bound = distances.iter().copied().fold(0.0_f64, f64::max) + allowance;
        if !bound.is_finite()
            || distances[0] + allowance > self.budget
            || distances[distances.len() - 1] + allowance > self.budget
        {
            return None;
        }
        if bound <= self.budget && self.monotone(&controls_vec, allowance, angular) {
            return Some(bound);
        }
        if depth == 0 {
            return None;
        }
        let (left, right) = subdivide_homogeneous_bezier(controls)?;
        Some(
            self.certify(&left, allowance, depth - 1, remaining, angular)?
                .max(self.certify(&right, allowance, depth - 1, remaining, angular)?),
        )
    }

    fn certify_curve(&self, curve: &BsplineCurve<Vector4>) -> Option<FacetCurveCertificate> {
        let controls = curve.control_points();
        let min_weight = controls.iter().map(|p| p.w).fold(f64::INFINITY, f64::min);
        let max_weight = controls.iter().map(|p| p.w).fold(0.0_f64, f64::max);
        let scale = controls
            .iter()
            .map(|&p| rational_point(p))
            .chain(self.source)
            .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
            .fold(1.0_f64, f64::max);
        // Same conditioning convention as the existing segment locus
        // certificate. This is a numerical margin, not directed rounding.
        let degree = crate::usize_to_f64(curve.degree());
        let allowance =
            256.0 * f64::EPSILON * (degree + 1.0).powi(2) * scale * (max_weight / min_weight);
        if !allowance.is_finite() || allowance >= self.budget {
            return None;
        }
        let spans = proof_bezier_spans(curve)?;
        let point = |p: [f64; 4]| rational_point(Vector4::new(p[0], p[1], p[2], p[3]));
        let front = point(*spans.first()?.controls.first()?);
        let back = point(*spans.last()?.controls.last()?);
        let actual_rays = [front - self.apex, back - self.apex];
        let seam_gaps = std::array::from_fn(|i| (actual_rays[i] - self.rays[i]).magnitude());
        let project = |ray: Vector3| ray - self.normal * self.normal.dot(ray);
        let [a, b] = actual_rays.map(project);
        let angular_error = (64.0 * f64::EPSILON * a.magnitude()).mul_add(
            b.magnitude(),
            allowance * (a.magnitude() + b.magnitude() + allowance),
        );
        if cross(self.normal, a, b) <= angular_error {
            return None;
        }
        let shared_rays = self.rays.map(project);
        if shared_rays.iter().any(|ray| ray.magnitude() <= allowance) {
            return None;
        }
        // A positive half-plane, plus strict angular progress, prevents
        // the curve from wrapping through the apex or leaving its actual
        // endpoint wedge. It does not depend on nearly equal seam rays.
        let half_plane = a.normalize() + b.normalize();
        if half_plane.magnitude() <= 64.0 * f64::EPSILON {
            return None;
        }
        let seam_deviation =
            2.0_f64.mul_add(allowance, seam_gaps.into_iter().fold(0.0_f64, f64::max));
        if !seam_deviation.is_finite() || seam_deviation > crate::KERNEL_EPSILON.min(self.budget) {
            return None;
        }
        let angular = AngularWedge {
            half_plane: half_plane.normalize(),
            shared_rays,
            seam_gaps,
        };
        let endpoint_bound = (front - self.source[0])
            .magnitude()
            .max((back - self.source[1]).magnitude())
            + allowance;
        if endpoint_bound > self.budget {
            return None;
        }
        // A positive-width knot interval may contain a discontinuity if
        // interior multiplicity is degree+1. The fan is a continuous map,
        // so reject unresolved jumps instead of treating them as chords.
        for pair in spans.windows(2) {
            let left = point(*pair[0].controls.last()?);
            let right = point(*pair[1].controls.first()?);
            if (left - right).magnitude() > allowance {
                return None;
            }
        }
        let mut remaining = 32768;
        let hull = spans.iter().try_fold(0.0_f64, |bound, span| {
            Some(bound.max(self.certify(
                &span.controls,
                allowance,
                16,
                &mut remaining,
                &angular,
            )?))
        })?;
        // Include the maximum numerically unresolved join gap for the
        // reverse (chord-to-curve) inclusion.
        let bound = hull.max(endpoint_bound) + allowance;
        (bound <= self.budget).then_some(FacetCurveCertificate {
            deviation: bound,
            seam_deviation,
        })
    }

    fn surface(&self, curve: &BsplineCurve<Vector4>) -> Option<Surface> {
        // Repeat the base weights in the collapsed row. This makes the
        // rational interpolation exactly (1-v)C(u)+v*apex, even for
        // nonuniform rational weights and multi-span boundary splines.
        let controls = curve
            .control_points()
            .iter()
            .map(|&p| vec![p, self.apex.to_vec().extend(1.0) * p.w])
            .collect();
        Some(Surface::NurbsSurface(NurbsSurface::new(
            BsplineSurface::new(
                (curve.knot_vector().clone(), KnotVector::bezier_knot(1)),
                controls,
            )
            .ok()?,
        )))
    }
}

/// Construct a boundary-constrained native approximation of a convex
/// source facet, from one oriented exact native edge per source side.
///
/// Each edge is certified against its source chord and its angular wedge.
/// The source and reconstructed fans have a common apex, initially the
/// source centroid. If a different interior apex is required, its
/// deviation from the reference source fan is charged to the budget. Continuity of
/// the curve's projection onto the chord, together with the endpoint
/// displacement bound, proves the reverse distance bound as well; convex
/// interpolation from the apex transfers both bounds to the whole fan.
/// The projected angular derivative prevents folds in each actual curve
/// wedge. The shared native wedge is additionally bounded by the measured
/// endpoint discrepancies, at most the kernel incidence tolerance; this
/// reports a seam uncertainty rather than claiming exact nonoverlap.
///
/// # Errors
/// Rejects invalid winding, degenerate or nonconvex projected facets,
/// incompatible native edge endpoints, unavailable exact rational nets,
/// and approximations that cannot meet the specified budget. No source
/// geometry or input edge identity is modified.
pub fn approximate_boundary_facet(
    source: &[Point3],
    edges: &[Edge],
    budget: f64,
) -> Result<BoundaryFacetApproximation, BoundaryFacetError> {
    use BoundaryFacetError::{InvalidBoundary, Uncertified, UnsupportedCurve};
    if source.len() < 3
        || source.len() != edges.len()
        || !budget.is_finite()
        || budget <= 0.0
        || !source.iter().copied().all(finite)
    {
        return Err(InvalidBoundary);
    }
    let apex = centroid(source);
    let normal = directed_pairs(source)
        .map(|(&a, &b)| (a - apex).cross(b - apex))
        .fold(Vector3::new(0.0, 0.0, 0.0), |sum, n| sum + n);
    let length = normal.magnitude();
    if length <= 0.0 || !length.is_finite() {
        return Err(InvalidBoundary);
    }
    let normal = normal / length;
    if !convex_fan(source, apex, normal)
        || !directed_pairs(edges).all(|(a, b)| a.back() == b.front())
    {
        return Err(InvalidBoundary);
    }
    let native = edges.iter().map(|e| e.front().point()).collect::<Vec<_>>();
    if !native.iter().copied().all(finite)
        || !source
            .iter()
            .zip(&native)
            .all(|(&a, &b)| (a - b).magnitude() <= budget)
    {
        return Err(InvalidBoundary);
    }
    let curves = edges
        .iter()
        .map(|edge| {
            if !edge.is_geometric_consistent() {
                return Err(InvalidBoundary);
            }
            let curve = edge
                .oriented_curve()
                .try_into_homogeneous_bspline_curve()
                .ok_or(UnsupportedCurve)?;
            if !positive_weights(curve.control_points()) {
                return Err(UnsupportedCurve);
            }
            Ok(curve)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let fan = FacetFan {
        source,
        edges,
        native: &native,
        curves: &curves,
        normal,
        budget,
    };
    match fan.construct(apex, 0.0) {
        Ok(result) => Ok(result),
        Err(error @ (InvalidBoundary | Uncertified)) => fan.recover(apex, error),
        Err(error) => Err(error),
    }
}

struct FacetFan<'a> {
    source: &'a [Point3],
    edges: &'a [Edge],
    native: &'a [Point3],
    curves: &'a [BsplineCurve<Vector4>],
    normal: Vector3,
    budget: f64,
}

impl FacetFan<'_> {
    fn construct(
        &self,
        apex: Point3,
        source_penalty: f64,
    ) -> Result<BoundaryFacetApproximation, BoundaryFacetError> {
        use BoundaryFacetError::{InvalidBoundary, NativeFace, Uncertified};
        if !source_penalty.is_finite()
            || source_penalty >= self.budget
            || !convex_fan(self.source, apex, self.normal)
            || !convex_fan(self.native, apex, self.normal)
        {
            return Err(InvalidBoundary);
        }
        let apex_vertex = builder::vertex(apex);
        let radial = self
            .edges
            .iter()
            .map(|edge| builder::line(edge.front(), &apex_vertex))
            .collect::<Vec<_>>();
        let mut faces = Vec::with_capacity(self.edges.len());
        let mut upper = 0.0_f64;
        let mut seam_upper = 0.0_f64;
        for (i, edge) in self.edges.iter().enumerate() {
            let next = (i + 1) % self.edges.len();
            let curve = &self.curves[i];
            let wedge = Wedge {
                apex,
                normal: self.normal,
                source: [self.source[i], self.source[next]],
                rays: [self.native[i] - apex, self.native[next] - apex],
                budget: self.budget - source_penalty,
            };
            let certificate = wedge.certify_curve(curve).ok_or(Uncertified)?;
            upper = upper
                .max(certificate.deviation)
                .max((self.native[i] - self.source[i]).magnitude());
            seam_upper = seam_upper.max(certificate.seam_deviation);
            let surface = wedge.surface(curve).ok_or(NativeFace)?;
            let face = Face::new(
                vec![Wire::from(vec![
                    edge.clone(),
                    radial[next].clone(),
                    radial[i].inverse(),
                ])],
                surface,
            )
            .map_err(|_| NativeFace)?;
            if !face.is_geometric_consistent() {
                return Err(NativeFace);
            }
            faces.push(face);
        }
        if upper + source_penalty > self.budget {
            return Err(Uncertified);
        }
        Ok(BoundaryFacetApproximation {
            faces,
            deviation_upper_bound: upper + source_penalty,
            seam_deviation_upper_bound: seam_upper,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Curve, ParametricSurface, Shell};

    fn source() -> [Point3; 3] {
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ]
    }

    fn edges(curvature: f64) -> Result<Vec<Edge>, Box<dyn std::error::Error>> {
        let [a, b, c] = source().map(builder::vertex);
        let curve = Curve::NurbsCurve(crate::NurbsCurve::new(BsplineCurve::new(
            KnotVector::bezier_knot(2),
            vec![
                a.point().to_vec().extend(1.0),
                Vector4::new(0.7, 0.0, curvature * 0.7, 0.7),
                b.point().to_vec().extend(1.0),
            ],
        )?));
        Ok(vec![
            Edge::new(&a, &b, curve)?,
            builder::line(&b, &c),
            builder::line(&c, &a),
        ])
    }

    #[test]
    fn curved_boundary_can_require_a_noncentroid_apex() -> Result<(), Box<dyn std::error::Error>> {
        let points = source();
        let [a, b, c] = points.map(builder::vertex);
        let edge = Edge::new(
            &a,
            &b,
            Curve::NurbsCurve(crate::NurbsCurve::new(BsplineCurve::new(
                KnotVector::bezier_knot(2),
                vec![
                    a.point().to_vec().extend(1.0),
                    Vector4::new(0.2, 0.5, 0.0, 1.0),
                    b.point().to_vec().extend(1.0),
                ],
            )?)),
        )?;
        let edges = [edge, builder::line(&b, &c), builder::line(&c, &a)];
        let curves = edges
            .iter()
            .map(|edge| {
                edge.oriented_curve()
                    .try_into_homogeneous_bspline_curve()
                    .ok_or("curve")
            })
            .collect::<Result<Vec<_>, _>>()?;
        let fan = FacetFan {
            source: &points,
            edges: &edges,
            native: &points,
            curves: &curves,
            normal: Vector3::new(0.0, 0.0, 1.0),
            budget: 0.3,
        };
        assert!(matches!(
            fan.construct(centroid(&points), 0.0),
            Err(BoundaryFacetError::Uncertified)
        ));
        let approximation = approximate_boundary_facet(&points, &edges, 0.3)?;
        assert!(approximation.deviation_upper_bound <= 0.3);
        assert!(
            approximation
                .faces
                .iter()
                .all(Face::is_geometric_consistent)
        );
        assert_eq!(approximation.faces[0].boundaries()[0][0], edges[0]);
        assert!(approximate_boundary_facet(&points, &edges, 0.1).is_err());
        Ok(())
    }

    #[test]
    fn near_tangent_endpoint_discrepancy_is_bounded_not_snapped()
    -> Result<(), Box<dyn std::error::Error>> {
        let [a, b, c] = source().map(builder::vertex);
        let offset = 5.0e-8;
        let curve = Curve::NurbsCurve(crate::NurbsCurve::new(BsplineCurve::new(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(-offset, 0.0, 0.0, 1.0),
                Vector4::new(1.0, 0.0, 0.05, 1.0),
                b.point().to_vec().extend(1.0),
            ],
        )?));
        let edges = [
            Edge::new(&a, &b, curve)?,
            builder::line(&b, &c),
            builder::line(&c, &a),
        ];
        let approximation = approximate_boundary_facet(&source(), &edges, 0.1)?;
        assert!(approximation.seam_deviation_upper_bound >= offset);
        assert!(approximation.seam_deviation_upper_bound < 1.0e-7);
        let wire = approximation.faces[0]
            .boundaries()
            .into_iter()
            .next()
            .ok_or("boundary")?;
        assert_eq!(wire[0], edges[0]);
        assert!(
            approximation
                .faces
                .iter()
                .all(Face::is_geometric_consistent)
        );
        Ok(())
    }

    #[test]
    fn rational_curved_edge_is_shared_and_whole_facet_is_bounded()
    -> Result<(), Box<dyn std::error::Error>> {
        let edges = edges(0.2)?;
        let approximation = approximate_boundary_facet(&source(), &edges, 0.2)?;
        assert_eq!(approximation.faces.len(), 3);
        assert!(approximation.deviation_upper_bound <= 0.2);
        let shell = approximation.faces.iter().cloned().collect::<Shell>();
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Oriented
        );
        assert!(shell.is_geometric_consistent());
        for (i, face) in approximation.faces.iter().enumerate() {
            let wire = face.boundaries().into_iter().next().ok_or("boundary")?;
            assert_eq!(wire[0], edges[i]);
            for u in 0..=10 {
                for v in 0..=10 {
                    let p = face
                        .oriented_surface()
                        .evaluate(f64::from(u) / 10.0, f64::from(v) / 10.0);
                    assert!(p.z.abs() <= approximation.deviation_upper_bound + 1.0e-12);
                }
            }
        }
        assert!(approximate_boundary_facet(&source(), &edges, 0.001).is_err());
        Ok(())
    }
}
