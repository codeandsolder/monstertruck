//! Marching intersection-curve backend, descended from the published 0.3.2 kernel.
//!
//! It still extracts mesh-vs-mesh interference segments and re-samples them onto
//! both carrier surfaces, but endpoint stitching is now fail-closed around branch
//! ambiguity and uses explicit sampled validation for short linear gaps. The old
//! parameter-space wrapper remains intentionally absent; the loops store consumes
//! only the cleaned 3-D [`IntersectionCurve`] leader.

use monstertruck_core::{
    cgmath64::{EuclideanSpace, InnerSpace, MetricSpace, Point3},
    tolerance::{DEFAULT_NUMERICAL_TOLERANCES, KERNEL_EPSILON},
};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use num_traits::ToPrimitive;

use crate::transversal::polyline_construction::{
    construct_polylines, enumerate_endpoint_joins, merge_endpoint_join,
    remove_redundant_interference_strands, stitch_nearby_polylines,
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

/// Track the exact UV branch along adjacent samples of one short SSI join.
/// Each candidate is still independently checked on the carrier; the hint
/// skips repeated global NURBS presearch when consecutive samples belong to
/// the same smooth trim. Fall back to the original unhinted projection whenever
/// local continuation cannot establish point-on-surface incidence.
fn point_is_on_surface<S>(
    surface: &S,
    point: Point3,
    tolerance: f64,
    hint: &mut Option<(f64, f64)>,
) -> bool
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let within_tolerance = |parameter: (f64, f64)| {
        surface.evaluate(parameter.0, parameter.1).distance(point) <= tolerance
    };
    let local = hint.and_then(|previous| {
        surface
            .search_parameter(point, Some(previous), 100)
            .filter(|parameter| within_tolerance(*parameter))
            .or_else(|| {
                surface
                    .search_nearest_parameter(point, Some(previous), 100)
                    .filter(|parameter| within_tolerance(*parameter))
            })
    });
    let parameter = local.or_else(|| {
        surface
            .search_parameter(point, None, 100)
            .filter(|parameter| within_tolerance(*parameter))
            .or_else(|| {
                surface
                    .search_nearest_parameter(point, None, 100)
                    .filter(|parameter| within_tolerance(*parameter))
            })
    });
    if let Some(parameter) = parameter {
        *hint = Some(parameter);
        true
    } else {
        false
    }
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
        .max(2.0);
    if !intervals.is_finite() || intervals > 4096.0 {
        return false;
    }
    let samples: Vec<_> = (0_u32..=4096)
        .map(f64::from)
        .take_while(|&index| index <= intervals)
        .map(|index| a + bridge * (index / intervals))
        .collect();
    // The carrier polygons are tessellations, not proof of exact incidence,
    // but their trim bounds can cheaply reject impossible candidates *before*
    // the costly global NURBS inverse search. This only reorders two required
    // predicates, preserving the exact admission criteria.
    let polygon_tolerance = 2.0 * validation.tolerance;
    if !validation
        .polygon0
        .neighborhood_include(&samples, polygon_tolerance)
        || !validation
            .polygon1
            .neighborhood_include(&samples, polygon_tolerance)
    {
        return false;
    }
    let mut hint0 = None;
    let mut hint1 = None;
    samples.iter().copied().all(|point| {
        point_is_on_surface(validation.surface0, point, validation.tolerance, &mut hint0)
            && point_is_on_surface(validation.surface1, point, validation.tolerance, &mut hint1)
    })
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
            let join_endpoints = [join.first, join.second];
            let shares_endpoint =
                join_endpoints.contains(&other.first) || join_endpoints.contains(&other.second);
            let distance_delta = (other.distance2.sqrt() - join_distance).abs();
            shares_endpoint && distance_delta <= tolerance
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
    // The original rescanning algorithm discarded the exact UV solution at
    // every sampled SSI vertex. On NURBS supports this makes the *next* point
    // redo two global 50x50 parameter preselections, before a 4-variable solve.
    // Continue the same intersection branch using the previously certified UVs
    // as Newton hints. Fall back to the existing global solver if local
    // continuation cannot prove exact incidence near the sampled leader.
    let mut prior = None::<((f64, f64), (f64, f64))>;
    let mut polyline = PolylineCurve(Vec::with_capacity(len));
    let support_tol =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_surface_projection_factor * KERNEL_EPSILON;
    let leader_tol = DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor;
    for index in 0..len {
        if index == len - 1 && raw[0].near(&raw[index]) {
            let first = polyline[0];
            polyline.push(first);
            break;
        }
        let parameter = index.to_f64()?;
        let point = raw[index];
        let tangent = ic.leader().derivative(parameter);
        let continuation = prior
            .filter(|_| tangent.magnitude2() > KERNEL_EPSILON * KERNEL_EPSILON)
            .and_then(|(uv0, uv1)| {
                ic.search_intersection_in_plane(point, tangent, Some(uv0), Some(uv1), 100)
            })
            .filter(|(candidate, uv0, uv1)| {
                // Newton convergence alone is not a correctness certificate.
                // Both carriers must agree on the SAME 3-D point and the
                // candidate must remain on the mesh-identified local branch.
                ic.surface0()
                    .evaluate(uv0.x, uv0.y)
                    .distance(ic.surface1().evaluate(uv1.x, uv1.y))
                    <= support_tol
                    && candidate.distance(point) <= leader_tol
            });
        let (q, uv0, uv1) = continuation
            .or_else(|| ic.search_triple(parameter, 100))
            .or_else(|| recover_triple(point))?;
        prior = Some(((uv0.x, uv0.y), (uv1.x, uv1.y)));
        polyline.push(q);
    }
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
    let stitch_tolerance = 4.0 * tol.max(KERNEL_EPSILON);
    let polylines =
        remove_redundant_interference_strands(construct_polylines(&interferences), KERNEL_EPSILON);
    let polylines = stitch_nearby_polylines(polylines, stitch_tolerance)?;
    // This pass repairs numerical fragmentation introduced by tessellation.
    // Its admission envelope is an explicit kernel-precision bound rather than
    // a threshold derived from the user-facing modeling-error budget.
    let validation_tolerance =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_fragment_validation_factor * KERNEL_EPSILON;
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
        .expect("fixture mesh is valid")
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
        .expect("fixture mesh is valid")
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
        .expect("fixture mesh is valid")
    }

    #[test]
    fn intersection_leader_continuation_certifies_dense_transverse_planes() {
        let leader = PolylineCurve(
            (0..=100)
                .map(|i| {
                    let x = -0.5 + f64::from(i) / 100.0;
                    Point3::new(x, 0.0, 0.0)
                })
                .collect(),
        );
        let recovered = build_intersection_curve(Plane::xy(), Plane::zx(), leader)
            .expect("transverse planes share an exact line");
        assert_eq!(recovered.leader().len(), 101);
        for (i, point) in recovered.leader().iter().enumerate() {
            assert!((point.x - (-0.5 + i as f64 / 100.0)).abs() < 1.0e-6);
            assert!(point.y.abs() < 1.0e-6);
            assert!(point.z.abs() < 1.0e-6);
        }
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
