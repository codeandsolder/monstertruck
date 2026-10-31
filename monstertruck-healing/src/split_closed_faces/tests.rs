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
