//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::split_closed_faces;
use monstertruck_geometry::prelude::{
    BoundedCurve, Cut, EuclideanSpace, Line, Matrix4, ParameterBoundary2D, ParameterCurve,
    ParameterDivision1D, ParametricCurve, Point2, Point3, Processor, RevolutionSurface,
    SearchNearestParameterD1, SearchParameter, Transformed, TrimmedCurve, UnitCircle, Vector3,
};
use monstertruck_topology::{
    Shell, StableId,
    compress::{CompressedEdge, CompressedEdgeIndex, CompressedFace as Face, CompressedShell},
};
use std::f64::consts::PI;

type Surface = RevolutionSurface<Line<Point3>>;

#[derive(
    Clone, Debug, ParametricCurve, BoundedCurve, ParameterDivision1D, Cut, SearchNearestParameterD1,
)]
enum Curve {
    Line(Line<Point3>),
    Arc(TrimmedCurve<Processor<UnitCircle<Point3>, Matrix4>>),
    ParameterCurve(ParameterCurve<Line<Point2>, Surface>),
}

impl From<ParameterCurve<Line<Point2>, Surface>> for Curve {
    fn from(value: ParameterCurve<Line<Point2>, Surface>) -> Self { Self::ParameterCurve(value) }
}

impl ParameterBoundary2D<Surface> for Curve {
    fn parameter_boundary_2d(&self, surface: &Surface, tolerance: f64) -> Option<Vec<Point2>> {
        match self {
            Self::Line(curve) => curve.parameter_boundary_2d(surface, tolerance),
            Self::Arc(curve) => curve.parameter_boundary_2d(surface, tolerance),
            Self::ParameterCurve(curve) => curve.parameter_boundary_2d(surface, tolerance),
        }
    }
}

fn sp(surface: &Surface, point: Point3, hint: Option<(f64, f64)>) -> Option<(f64, f64)> {
    surface.search_parameter(point, hint, 10)
}

fn simple_cylinder_shell() -> CompressedShell<Point3, Curve, Surface> {
    let vertices = vec![
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(-1.0, 0.0, 1.0),
    ];
    let translate = Matrix4::from_translation(Vector3::unit_z());
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: Curve::Arc(TrimmedCurve::new(
                Processor::new(UnitCircle::new()),
                (0.0, PI),
            )),
        },
        CompressedEdge {
            vertices: (1, 0),
            curve: Curve::Arc(TrimmedCurve::new(
                Processor::new(UnitCircle::new()),
                (PI, 2.0 * PI),
            )),
        },
        CompressedEdge {
            vertices: (0, 2),
            curve: Curve::Line(Line(vertices[0], vertices[2])),
        },
        CompressedEdge {
            vertices: (2, 3),
            curve: Curve::Arc(TrimmedCurve::new(
                Processor::new(UnitCircle::new()).transformed(translate),
                (0.0, PI),
            )),
        },
        CompressedEdge {
            vertices: (3, 2),
            curve: Curve::Arc(TrimmedCurve::new(
                Processor::new(UnitCircle::new()).transformed(translate),
                (PI, 2.0 * PI),
            )),
        },
    ];
    let surface = RevolutionSurface::by_revolution(
        Line(vertices[2], vertices[0]),
        Point3::origin(),
        Vector3::unit_z(),
    );
    let faces = vec![Face {
        surface,
        boundaries: vec![vec![
            CompressedEdgeIndex {
                index: 1,
                orientation: true,
            },
            CompressedEdgeIndex {
                index: 2,
                orientation: true,
            },
            CompressedEdgeIndex {
                index: 4,
                orientation: false,
            },
            CompressedEdgeIndex {
                index: 3,
                orientation: false,
            },
            CompressedEdgeIndex {
                index: 2,
                orientation: false,
            },
            CompressedEdgeIndex {
                index: 0,
                orientation: true,
            },
        ]],
        orientation: true,
    }];
    CompressedShell {
        vertices,
        edges,
        faces,
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    }
}

#[test]
fn split_closed_faces_reports_change_only_on_first_normalizing_pass() {
    let mut shell = simple_cylinder_shell();
    let original_vertex_count = shell.vertices.len();
    let original_edge_count = shell.edges.len();
    let original_face_count = shell.faces.len();
    shell.vertex_stable_ids = Some(
        (100_u64..)
            .take(original_vertex_count)
            .map(StableId::new)
            .collect(),
    );
    shell.edge_stable_ids = Some(
        (200_u64..)
            .take(original_edge_count)
            .map(StableId::new)
            .collect(),
    );
    shell.face_stable_ids = Some(
        (300_u64..)
            .take(original_face_count)
            .map(StableId::new)
            .collect(),
    );

    assert!(split_closed_faces(&mut shell, 0.01, sp));
    assert!(Shell::extract(shell.clone()).is_ok());
    let vertex_ids = shell.vertex_stable_ids.as_ref().expect("vertex stable ids");
    let edge_ids = shell.edge_stable_ids.as_ref().expect("edge stable ids");
    let face_ids = shell.face_stable_ids.as_ref().expect("face stable ids");
    assert_eq!(vertex_ids.len(), shell.vertices.len());
    assert_eq!(edge_ids.len(), shell.edges.len());
    assert_eq!(face_ids.len(), shell.faces.len());
    assert_eq!(
        &vertex_ids[..original_vertex_count],
        &(100_u64..)
            .take(original_vertex_count)
            .map(StableId::new)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        &edge_ids[..original_edge_count],
        &(200_u64..)
            .take(original_edge_count)
            .map(StableId::new)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        &face_ids[..original_face_count],
        &(300_u64..)
            .take(original_face_count)
            .map(StableId::new)
            .collect::<Vec<_>>()
    );
    assert!(
        vertex_ids[original_vertex_count..]
            .iter()
            .all(|id| *id == StableId::UNASSIGNED)
    );
    assert!(
        edge_ids[original_edge_count..]
            .iter()
            .all(|id| *id == StableId::UNASSIGNED)
    );
    assert!(
        face_ids[original_face_count..]
            .iter()
            .all(|id| *id == StableId::UNASSIGNED)
    );

    let stabilized = shell.clone();
    assert!(!split_closed_faces(&mut shell, 0.01, sp));
    assert_eq!(shell.vertices, stabilized.vertices);
    assert_eq!(shell.edges.len(), stabilized.edges.len());
    assert_eq!(shell.faces.len(), stabilized.faces.len());
}

#[test]
fn healing_boundary_closure_uses_each_axis_native_scale_and_3d_evidence() {
    use super::{BoundaryType, PolylineCurve, boundary_type};
    use monstertruck_geometry::prelude::{BsplineSurface, KnotVector};

    let square_surface = |u_span: f64, v_span: f64, physical_u_span: f64| {
        BsplineSurface::new_unchecked(
            (
                KnotVector::from(vec![0.0, 0.0, u_span, u_span]),
                KnotVector::from(vec![0.0, 0.0, v_span, v_span]),
            ),
            vec![
                vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
                vec![
                    Point3::new(physical_u_span, 0.0, 0.0),
                    Point3::new(physical_u_span, 1.0, 0.0),
                ],
            ],
        )
    };
    for (u_span, v_span) in [(1e-9, 1e9), (1e9, 1e-9), (1.0, 1.0)] {
        let surface = square_surface(u_span, v_span, 1.0);
        let close = Point2::new(0.5e-6 * u_span, 0.5e-6 * v_span);
        let poly = PolylineCurve(vec![
            Point2::new(0.0, 0.0),
            Point2::new(u_span, 0.0),
            Point2::new(u_span, v_span),
            Point2::new(0.0, v_span),
            close,
        ]);
        assert!(
            matches!(boundary_type(&poly, &surface, 1e-4), BoundaryType::Positive),
            "physical unit-square should close regardless of native U/V scale {u_span:e}/{v_span:e}"
        );
        let reversed = PolylineCurve(vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.0, v_span),
            Point2::new(u_span, v_span),
            Point2::new(u_span, 0.0),
            close,
        ]);
        assert!(matches!(
            boundary_type(&reversed, &surface, 1e-4),
            BoundaryType::Negative
        ));
        let interrupted = PolylineCurve(vec![
            Point2::new(0.0, 0.0),
            Point2::new(u_span, 0.0),
            Point2::new(u_span, v_span),
            Point2::new(0.0, v_span),
            Point2::new(1e-4 * u_span, 0.0),
        ]);
        assert!(
            matches!(
                boundary_type(&interrupted, &surface, 1e-4),
                BoundaryType::NotClosed
            ),
            "a 0.01% native UV gap should never be closed by an absolute 1e-3 squared-UV cutoff"
        );
    }
    let large_physical = square_surface(1.0, 1.0, 1.0e9);
    let uv_small_but_physically_open = PolylineCurve(vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
        Point2::new(0.5e-6, 0.0),
    ]);
    assert!(
        matches!(
            boundary_type(&uv_small_but_physically_open, &large_physical, 1e-4),
            BoundaryType::NotClosed
        ),
        "small normalized UV residual alone is not a physically closed wire"
    );
}

#[test]
fn healing_boundary_closure_rejects_nonfinite_and_short_contours() {
    use super::{BoundaryType, PolylineCurve, boundary_type};
    use monstertruck_geometry::prelude::Plane;

    let plane = Plane::xy();
    let short = PolylineCurve(vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.0, 0.0),
    ]);
    assert!(matches!(
        boundary_type(&short, &plane, 1e-6),
        BoundaryType::NotClosed
    ));
    let invalid = PolylineCurve(vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, f64::NAN),
        Point2::new(0.0, 0.0),
    ]);
    assert!(matches!(
        boundary_type(&invalid, &plane, 1e-6),
        BoundaryType::NotClosed
    ));
    let valid = PolylineCurve(vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.0, 1.0),
        Point2::new(0.0, 0.0),
    ]);
    assert!(matches!(
        boundary_type(&valid, &plane, f64::INFINITY),
        BoundaryType::NotClosed
    ));
    assert!(matches!(
        boundary_type(&valid, &plane, -1.0),
        BoundaryType::NotClosed
    ));
    let degenerate = PolylineCurve(vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(0.0, 0.0),
    ]);
    assert!(matches!(
        boundary_type(&degenerate, &plane, 1.0e-6),
        BoundaryType::NotClosed
    ));
}

#[test]
fn normalized_concave_divisor_is_invariant_to_independent_uv_affine_scales() {
    use super::{concave_wire_divisor, normalized_polygon, polygon_contains_point};
    let vertices = [
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 1.0),
        (2.0, 1.0),
        (2.0, 3.0),
        (0.0, 3.0),
    ];
    let mut expected_diagonal = None;
    for (u_span, v_span) in [
        (1.0e-9, 1.0e9),
        (1.0e9, 1.0e-9),
        (1.0e-6, 1.0e-6),
        (1.0, 1.0),
        (1.0e9, 1.0e9),
    ] {
        let polygon = vertices
            .into_iter()
            .map(|(u, v)| Point2::new(u / 4.0 * u_span, v / 3.0 * v_span))
            .collect::<Vec<_>>();
        let indexed = polygon
            .iter()
            .copied()
            .enumerate()
            .map(|(i, point)| (i + 10, point))
            .collect::<Vec<_>>();
        let result = concave_wire_divisor(&indexed);
        assert!(
            result.is_some(),
            "lost concavity at native U/V {u_span:e}/{v_span:e}"
        );
        if let Some(expected) = expected_diagonal {
            assert_eq!(
                result,
                Some(expected),
                "different divisor for identical polygon"
            );
        } else {
            expected_diagonal = result;
        }

        let normalized = normalized_polygon(&polygon).expect("positive finite chart");
        assert!(polygon_contains_point(
            &normalized,
            Point2::new(0.125, 0.75)
        ));
        assert!(!polygon_contains_point(
            &normalized,
            Point2::new(0.75, 0.75)
        ));
        assert!(polygon_contains_point(
            &normalized,
            Point2::new(0.5, 1.0 / 3.0)
        ));
    }
}

#[test]
fn repeated_zero_length_segment_cannot_make_arbitrary_outside_point_inside() {
    use super::{polygon_contains_point, segment_contains_point};
    let vertex = Point2::new(0.0, 0.0);
    assert!(segment_contains_point(vertex, vertex, vertex));
    assert!(!segment_contains_point(
        vertex,
        vertex,
        Point2::new(0.75, 0.75)
    ));
    let polygon = vec![
        vertex,
        vertex,
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];
    assert!(!polygon_contains_point(&polygon, Point2::new(2.0, 2.0)));
    assert!(polygon_contains_point(&polygon, vertex));
    assert!(polygon_contains_point(&polygon, Point2::new(0.5, 0.5)));
}

#[test]
fn normalized_polygon_rejects_nonfinite_or_degenerate_charts() {
    use super::normalized_polygon;
    assert!(normalized_polygon(&[]).is_none());
    assert!(normalized_polygon(&[Point2::new(0.0, 0.0)]).is_none());
    assert!(normalized_polygon(&[Point2::new(0.0, 0.0), Point2::new(0.0, 2.0),]).is_none());
    assert!(normalized_polygon(&[Point2::new(0.0, 0.0), Point2::new(f64::NAN, 1.0),]).is_none());
    assert!(
        normalized_polygon(&[Point2::new(0.0, 0.0), Point2::new(f64::INFINITY, 1.0),]).is_none()
    );
}

#[test]
fn diagonal_strict_intersection_is_not_an_endpoint_touch() {
    use super::segments_intersect_strict;
    let a = Point2::new(0.0, 0.0);
    let b = Point2::new(1.0, 1.0);
    assert!(segments_intersect_strict(
        a,
        b,
        Point2::new(0.0, 1.0),
        Point2::new(1.0, 0.0)
    ));
    assert!(!segments_intersect_strict(a, b, b, Point2::new(1.0, 0.0)));
    assert!(!segments_intersect_strict(
        a,
        b,
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 1.0)
    ));
}

#[test]
fn seam_section_is_parameter_and_physical_scale_aware() {
    use super::{PolylineCurve, same_nonperiodic_section};
    use monstertruck_geometry::prelude::{BsplineSurface, KnotVector};

    let plane = |u_span: f64, v_span: f64, world_u_span: f64| {
        BsplineSurface::new_unchecked(
            (
                KnotVector::from(vec![0.0, 0.0, u_span, u_span]),
                KnotVector::from(vec![0.0, 0.0, v_span, v_span]),
            ),
            vec![
                vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
                vec![
                    Point3::new(world_u_span, 0.0, 0.0),
                    Point3::new(world_u_span, 1.0, 0.0),
                ],
            ],
        )
    };
    for (u_span, v_span) in [(1e-9, 1e9), (1e9, 1e-9), (1.0, 1.0)] {
        let surface = plane(u_span, v_span, 1.0);
        let p = Point2::new(0.0, 0.0);
        let q = Point2::new(u_span * 0.5e-6, v_span);
        let wire = vec![PolylineCurve(vec![p, Point2::new(u_span, v_span / 2.0), q])];
        assert!(
            same_nonperiodic_section(&surface, &wire, p, q, 0, 1e-4),
            "native U near-constant section failed with spans {u_span:e}/{v_span:e}"
        );
        assert!(!same_nonperiodic_section(&surface, &wire, p, q, 1, 1e-4));
        let far_q = Point2::new(u_span * 0.01, v_span);
        assert!(!same_nonperiodic_section(
            &surface, &wire, p, far_q, 0, 1e-4
        ));

        let v_small = Point2::new(u_span, v_span * 0.5e-6);
        assert!(same_nonperiodic_section(
            &surface,
            &[PolylineCurve(vec![
                p,
                Point2::new(u_span / 2.0, v_span),
                v_small
            ])],
            p,
            v_small,
            1,
            1e-4,
        ));
    }
    let huge_world = plane(1.0, 1.0, 1e9);
    let p = Point2::new(0.0, 0.0);
    let q = Point2::new(0.5e-6, 1.0);
    let wire = vec![PolylineCurve(vec![p, Point2::new(1.0, 0.5), q])];
    assert!(
        !same_nonperiodic_section(&huge_world, &wire, p, q, 0, 1e-4),
        "near in UV is not near in 3D on a huge surface"
    );
    assert!(!same_nonperiodic_section(&huge_world, &wire, p, q, 2, 1e-4));
    assert!(!same_nonperiodic_section(&huge_world, &[], p, q, 0, 1e-4));
    assert!(!same_nonperiodic_section(
        &huge_world,
        &wire,
        p,
        q,
        0,
        f64::NAN
    ));
}

#[test]
fn closest_native_period_representation_has_no_arbitrary_turn_limit() {
    use super::get_mindiff;
    for period in [1.0e-9, 1.0, 1.0e9] {
        let displaced = 10_001.25 * period;
        let matched = get_mindiff(displaced, 0.0, period).expect("valid positive period");
        assert!((matched - 0.25 * period).abs() <= 1.0e-8 * period);
        let opposite =
            get_mindiff(-10_001.25 * period, 0.0, period).expect("valid positive period");
        assert!((opposite + 0.25 * period).abs() <= 1.0e-8 * period);
        assert_eq!(
            get_mindiff(0.25 * period, 0.25 * period, period),
            Some(0.25 * period)
        );
    }
    assert_eq!(get_mindiff(0.5, 0.0, 0.0), None);
    assert_eq!(get_mindiff(0.5, 0.0, f64::NAN), None);
    assert_eq!(get_mindiff(f64::INFINITY, 0.0, 1.0), None);
    assert_eq!(get_mindiff(0.5, f64::NAN, 1.0), None);
}
