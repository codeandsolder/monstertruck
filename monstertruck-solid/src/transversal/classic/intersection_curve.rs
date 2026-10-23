//! Marching intersection-curve backend, descended from the published 0.3.2 kernel.
//!
//! It still extracts mesh-vs-mesh interference segments and re-samples them onto
//! both carrier surfaces, but endpoint stitching is now fail-closed around branch
//! ambiguity and uses explicit sampled validation for short linear gaps. The old
//! parameter-space wrapper remains intentionally absent; the loops store consumes
//! only the cleaned 3-D [`IntersectionCurve`] leader.

use monstertruck_core::{cgmath64::*, tolerance::TOLERANCE};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;

use crate::transversal::polyline_construction::{
    construct_polylines, enumerate_endpoint_joins, merge_endpoint_join, stitch_nearby_polylines,
};

type Polyline = PolylineCurve<Point3>;

struct FragmentValidation<'a, S> {
    surface0: &'a S,
    polygon0: &'a PolygonMesh,
    surface1: &'a S,
    polygon1: &'a PolygonMesh,
    tolerance: f64,
}

fn endpoint_and_inner(
    polyline: &Polyline,
    front: bool,
    minimum_span: f64,
) -> Option<(Point3, Point3)> {
    if polyline.len() < 2 {
        return None;
    }
    let endpoint = if front {
        polyline[0]
    } else {
        *polyline.last()?
    };
    let inner = if front {
        polyline
            .iter()
            .copied()
            .skip(1)
            .find(|point| endpoint.distance(*point) > minimum_span)?
    } else {
        polyline
            .iter()
            .copied()
            .rev()
            .skip(1)
            .find(|point| endpoint.distance(*point) > minimum_span)?
    };
    Some((endpoint, inner))
}

fn point_is_on_surface<S>(surface: &S, point: Point3, tolerance: f64) -> bool
where S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3> {
    let within_tolerance = |parameter: (f64, f64)| {
        surface.evaluate(parameter.0, parameter.1).distance(point) <= tolerance
    };
    if let Some(parameter) = surface.search_parameter(point, None, 100)
        && within_tolerance(parameter)
    {
        return true;
    }
    surface
        .search_nearest_parameter(point, None, 100)
        .is_some_and(within_tolerance)
}

fn validated_linear_fragment_join<S>(
    first: &Polyline,
    first_front: bool,
    second: &Polyline,
    second_front: bool,
    validation: &FragmentValidation<'_, S>,
) -> bool
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let Some((a, a_inner)) = endpoint_and_inner(first, first_front, validation.tolerance) else {
        return false;
    };
    let Some((b, b_inner)) = endpoint_and_inner(second, second_front, validation.tolerance) else {
        return false;
    };

    let bridge = b - a;
    let bridge_length = bridge.magnitude();

    let outward_a = a - a_inner;
    let outward_b = b - b_inner;
    let span_a = outward_a.magnitude();
    let span_b = outward_b.magnitude();
    if span_a <= validation.tolerance || span_b <= validation.tolerance {
        return false;
    }

    // The missing interval must continue outward from both fragments and remain
    // on their endpoint tangents to within the kernel validation tolerance.
    let direction_ok = outward_a.dot(bridge) > 0.0 && outward_b.dot(-bridge) > 0.0;
    let line_error_a = outward_a.cross(bridge).magnitude() / span_a;
    let line_error_b = outward_b.cross(bridge).magnitude() / span_b;
    let line_ok = line_error_a <= validation.tolerance && line_error_b <= validation.tolerance;
    if !(direction_ok && line_ok) {
        return false;
    }

    // Validate the entire short bridge, not just its endpoints. Dense samples must
    // stay on both exact carrier surfaces and inside both trimmed polygon faces.
    // Refuse pathological gaps rather than spending unbounded work validating them.
    let intervals = (bridge_length / (2.0 * validation.tolerance))
        .ceil()
        .max(2.0) as usize;
    if intervals > 4096 {
        return false;
    }
    let samples: Vec<_> = (0..=intervals)
        .map(|i| a + bridge * (i as f64 / intervals as f64))
        .collect();
    let surfaces_ok = samples.iter().copied().all(|point| {
        point_is_on_surface(validation.surface0, point, validation.tolerance)
            && point_is_on_surface(validation.surface1, point, validation.tolerance)
    });
    if !surfaces_ok {
        return false;
    }

    // The carrier polygons are tessellations, not the exact surfaces above.
    // Adjacent SSI paths already permit an 8x kernel-tolerance endpoint
    // canonicalization; use that same envelope here for trimmed-face membership
    // while keeping the exact surface/tangent validation at the tighter 4x bound.
    let polygon_tolerance = 2.0 * validation.tolerance;
    validation
        .polygon0
        .neighborhood_include(&samples, polygon_tolerance)
        && validation
            .polygon1
            .neighborhood_include(&samples, polygon_tolerance)
}

fn stitch_validated_linear_intersection_fragments<S>(
    mut polylines: Vec<Polyline>,
    surface0: &S,
    polygon0: &PolygonMesh,
    surface1: &S,
    polygon1: &PolygonMesh,
    tolerance: f64,
) -> Vec<Polyline>
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let validation = FragmentValidation {
        surface0,
        polygon0,
        surface1,
        polygon1,
        tolerance,
    };
    loop {
        let Some(mut joins) =
            enumerate_endpoint_joins(&polylines, |first, first_front, second, second_front| {
                validated_linear_fragment_join(
                    first,
                    first_front,
                    second,
                    second_front,
                    &validation,
                )
            })
        else {
            return polylines;
        };

        joins.sort_by(|left, right| {
            left.distance2
                .total_cmp(&right.distance2)
                .then_with(|| left.first.chain.cmp(&right.first.chain))
                .then_with(|| left.second.chain.cmp(&right.second.chain))
        });
        let Some(join) = joins.first().copied() else {
            return polylines;
        };

        // Several fragments may lie on one validated straight intersection, so an
        // endpoint can legitimately see more than one farther mate. Continue to
        // the closest one. A near-equal alternative sharing either endpoint is a
        // real branch ambiguity and remains fail-closed.
        let join_distance = join.distance2.sqrt();
        let ambiguous = joins.iter().skip(1).any(|other| {
            let shares_endpoint = other.first == join.first
                || other.second == join.first
                || other.first == join.second
                || other.second == join.second;
            shares_endpoint && (other.distance2.sqrt() - join_distance).abs() <= tolerance
        });
        if ambiguous {
            return polylines;
        }

        if merge_endpoint_join(&mut polylines, join).is_none() {
            return polylines;
        }
    }
}

/// Re-sample a raw interference polyline onto both surfaces, returning a
/// cleaned polyline intersection curve. Mirrors 0.3.2
/// `IntersectionCurveWithParameters::try_new`, keeping only the 3D leader.
fn build_intersection_curve<S>(
    surface0: S,
    surface1: S,
    poly: Polyline,
) -> Option<IntersectionCurve<Polyline, S, S>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let ic = IntersectionCurve::new(surface0.clone(), surface1.clone(), poly);
    let raw = ic.leader().clone();
    let len = raw.len();
    if len < 2 {
        return None;
    }
    let recover_triple = |point: Point3| {
        ic.search_nearest_point(point, None, None, 100).or_else(|| {
            let p0 = ic
                .surface0()
                .search_parameter(point, None, 100)
                .or_else(|| ic.surface0().search_nearest_parameter(point, None, 100))?;
            let p1 = ic
                .surface1()
                .search_parameter(point, None, 100)
                .or_else(|| ic.surface1().search_nearest_parameter(point, None, 100))?;
            let q0 = ic.surface0().evaluate(p0.0, p0.1);
            let q1 = ic.surface1().evaluate(p1.0, p1.1);
            Some((q0.midpoint(q1), p0.into(), p1.into()))
        })
    };
    let mut polyline = PolylineCurve(Vec::new());
    for i in 0..len - 1 {
        let (q, _, _) = ic
            .search_triple(i as f64, 100)
            .or_else(|| recover_triple(raw[i]))?;
        polyline.push(q);
    }
    let q = if raw[0].near(&raw[len - 1]) {
        polyline[0]
    } else {
        ic.search_triple((len - 1) as f64, 100)
            .or_else(|| recover_triple(raw[len - 1]))?
            .0
    };
    polyline.push(q);
    Some(IntersectionCurve::new(surface0, surface1, polyline))
}

type IntersectionTuple<S> = (Polyline, IntersectionCurve<Polyline, S, S>);

/// Marching intersection curves between two trimmed surfaces, keyed by their
/// tessellated polygon meshes. Returns, per intersection branch, the raw
/// interference polyline and the cleaned polyline intersection curve.
pub(super) fn intersection_curves<S>(
    surface0: S,
    polygon0: &PolygonMesh,
    surface1: S,
    polygon1: &PolygonMesh,
    tol: f64,
) -> Option<Vec<IntersectionTuple<S>>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let interferences = polygon0.extract_interference(polygon1);
    let stitch_tolerance = 4.0 * tol.max(TOLERANCE);
    let polylines = stitch_nearby_polylines(construct_polylines(&interferences), stitch_tolerance)?;
    // This pass repairs numerical fragmentation introduced by tessellation.
    // Its admission envelope is a kernel-precision bound, not the caller's
    // potentially much coarser Boolean modeling tolerance.
    let validation_tolerance = 4.0 * TOLERANCE;
    let polylines = stitch_validated_linear_intersection_fragments(
        polylines,
        &surface0,
        polygon0,
        &surface1,
        polygon1,
        validation_tolerance,
    );
    polylines
        .into_iter()
        // Mesh interference can contain isolated point contacts. They do not
        // define a one-dimensional trim and therefore cannot divide either
        // face. Keep real intersection branches strict, but discard these
        // zero-dimensional artifacts instead of letting one make the whole
        // Boolean look like an SSI failure.
        .filter(|polyline| polyline.len() >= 2)
        .map(|polyline| {
            let curve =
                build_intersection_curve(surface0.clone(), surface1.clone(), polyline.clone())?;
            Some((polyline, curve))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xy_mesh(x0: f64, x1: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(x0, -1.0, 0.0),
                    Point3::new(x1, -1.0, 0.0),
                    Point3::new(x1, 1.0, 0.0),
                    Point3::new(x0, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
        )
    }

    fn disconnected_xy_mesh(gap: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(-1.0, -1.0, 0.0),
                    Point3::new(0.0, -1.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(-1.0, 1.0, 0.0),
                    Point3::new(gap, -1.0, 0.0),
                    Point3::new(2.0, -1.0, 0.0),
                    Point3::new(2.0, 1.0, 0.0),
                    Point3::new(gap, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]]),
        )
    }

    fn zx_mesh(x0: f64, x1: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(x0, 0.0, -1.0),
                    Point3::new(x1, 0.0, -1.0),
                    Point3::new(x1, 0.0, 1.0),
                    Point3::new(x0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
        )
    }

    #[test]
    fn validated_linear_stitch_skips_tiny_endpoint_segments() {
        let tolerance = 4.0e-6;
        let gap = 5.0e-5;
        let fragments = vec![
            PolylineCurve(vec![
                Point3::new(-0.5, 0.0, 0.0),
                Point3::new(-2.0e-6, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
            ]),
            PolylineCurve(vec![
                Point3::new(gap, 0.0, 0.0),
                Point3::new(gap + 2.0e-6, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ]),
        ];
        let stitched = stitch_validated_linear_intersection_fragments(
            fragments,
            &Plane::xy(),
            &xy_mesh(-1.0, 2.0),
            &Plane::zx(),
            &zx_mesh(-1.0, 2.0),
            tolerance,
        );
        assert_eq!(stitched.len(), 1);
        assert_eq!(stitched[0].front(), Point3::new(-0.5, 0.0, 0.0));
        assert_eq!(stitched[0].back(), Point3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn validated_linear_stitch_does_not_cross_real_trim_gap() {
        let tolerance = 4.0e-6;
        let gap = 5.0e-5;
        let fragments = vec![
            PolylineCurve(vec![
                Point3::new(-0.5, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
            ]),
            PolylineCurve(vec![Point3::new(gap, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]),
        ];
        let stitched = stitch_validated_linear_intersection_fragments(
            fragments,
            &Plane::xy(),
            &disconnected_xy_mesh(gap),
            &Plane::zx(),
            &zx_mesh(-1.0, 2.0),
            tolerance,
        );
        assert_eq!(stitched.len(), 2);
    }
}
