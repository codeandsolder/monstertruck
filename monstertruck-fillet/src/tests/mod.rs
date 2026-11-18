use itertools::Itertools;
use monstertruck_geometry::prelude::*;
use monstertruck_io::step::save::{StepCurve, StepFormat, StepLength, StepModels};
use monstertruck_meshing::prelude::*;

use super::geometry::*;
use super::types::*;

use monstertruck_traits::ParametricSurface;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

use super::{
    FilletError, FilletOptions, FilletProfile, RadiusSpec, fillet, fillet_along_wire, fillet_edges,
    fillet_edges_generic, fillet_with_side, topology,
};

mod accuracy;
mod edges;
mod profile;
mod radius;
mod surface;

// CSG/boolean-result fillet coverage lives in
// monstertruck-modeling/tests/fillet_test.rs (feature `fillet`): the
// modeling-curve conversions are feature-gated there, and enabling them in
// this crate's unit tests would pull in a second monstertruck-fillet whose
// fillet types cannot unify with the crate under test.

fn test_artifact_dir() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace_dir = manifest_dir.parent().unwrap_or_else(|| Path::new("."));
    let dir = workspace_dir.join("target").join("tests").join("fillet");
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn step_polyline(curve: &Curve) -> PolylineCurve<Point3> {
    let (t0, t1) = curve.range_tuple();
    let points = (0..=32)
        .map(|i| t0 + (t1 - t0) * (i as f64) / 32.0)
        .map(|t| curve.evaluate(t))
        .collect();
    PolylineCurve(points)
}

impl StepFormat for Curve {
    fn fmt(&self, idx: usize, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        StepFormat::fmt(&step_polyline(self), idx, f)
    }
}

impl StepLength for Curve {
    fn step_length(&self) -> usize { step_polyline(self).step_length() }
}

impl StepCurve for Curve {}

fn format_entity(entity: &(impl StepFormat + StepLength), idx: usize) -> (String, usize) {
    struct EntityFmt<'a, T>(&'a T, usize);

    impl<T: StepFormat> fmt::Display for EntityFmt<'_, T> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            StepFormat::fmt(self.0, self.1, f)
        }
    }

    let body = EntityFmt(entity, idx).to_string();
    (body, idx + entity.step_length())
}

fn format_step_file(name: &str, data_section: &str) -> String {
    format!(
        "ISO-10303-21;\n\
         HEADER;\n\
         FILE_DESCRIPTION(('Fillet Debug: {name}'), '2;1');\n\
         FILE_NAME('{name}.step', '', (''), (''), '', 'monstertruck', '');\n\
         FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n\
         ENDSEC;\n\
         DATA;\n\
         {data_section}\
         ENDSEC;\n\
         END-ISO-10303-21;\n"
    )
}

fn write_step(path: PathBuf, name: &str, data_section: String) {
    fs::write(path, format_step_file(name, &data_section)).unwrap();
}

fn dump_shell_step(name: &str, shell: &Shell) {
    let compressed = shell.compress();
    let mut models = StepModels::default();
    models.push_shell(&compressed);
    write_step(
        test_artifact_dir().join(format!("{name}.step")),
        name,
        models.to_string(),
    );
}

fn dump_surface_step<S>(name: &str, surface: &S)
where S: StepFormat + StepLength {
    let (entity, _) = format_entity(surface, 1);
    write_step(
        test_artifact_dir().join(format!("{name}.step")),
        name,
        entity,
    );
}

fn dump_surface_pair_step<S>(name: &str, surface0: &S, surface1: &S)
where S: StepFormat + StepLength {
    let (entity0, next) = format_entity(surface0, 1);
    let (entity1, _) = format_entity(surface1, next);
    write_step(
        test_artifact_dir().join(format!("{name}.step")),
        name,
        format!("{entity0}{entity1}"),
    );
}

/// Build a unit-height box with the first `FACE_COUNT` faces, in stable
/// top/front/right/back/left/bottom order, preserving shared edge identity.
fn build_box_with_faces<const FACE_COUNT: usize>() -> (Shell, [Edge; 12], Vec<Vertex>) {
    const XY_CORNERS: [(f64, f64); 4] = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let points: [Point3; 8] = std::array::from_fn(|index| {
        let (x, y) = XY_CORNERS[index % XY_CORNERS.len()];
        Point3::new(x, y, if index < 4 { 1.0 } else { 0.0 })
    });
    build_box_from_points::<FACE_COUNT>(points)
}

/// Build an indexed box topology from arbitrary vertex coordinates, reusing
/// the same 12 physical edges across all of its parametric planar faces.
fn build_box_from_points<const FACE_COUNT: usize>(
    points: [Point3; 8],
) -> (Shell, [Edge; 12], Vec<Vertex>) {
    const FACE_CORNERS: [[usize; 4]; 6] = [
        [0, 1, 2, 3], // 0: top
        [1, 0, 4, 5], // 1: front
        [2, 1, 5, 6], // 2: right
        [3, 2, 6, 7], // 3: back
        [0, 3, 7, 4], // 4: left
        [5, 4, 7, 6], // 5: bottom
    ];
    const EDGE_ENDPOINTS: [[usize; 2]; 12] = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [0, 4],
        [1, 5],
        [2, 6],
        [3, 7],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4],
    ];
    assert!(FACE_COUNT <= FACE_CORNERS.len());

    let vertices = Vertex::from_points(points);
    // Fixed indexing retains the edge-order contract used by every fillet test.
    let edges = EDGE_ENDPOINTS.map(|[front, back]| {
        let line = Line(points[front], points[back]);
        let curve = BsplineCurve::new_unchecked(KnotVector::bezier_knot(1), vec![line.0, line.1]);
        Edge::new_unchecked(
            &vertices[front],
            &vertices[back],
            NurbsCurve::from(curve).into(),
        )
    });
    let plane = |[i, j, k, l]: [usize; 4]| {
        let control_points = vec![vec![points[i], points[l]], vec![points[j], points[k]]];
        let knots = KnotVector::bezier_knot(1);
        let bspline = BsplineSurface::new_unchecked((knots.clone(), knots), control_points);
        let wire: Wire = [i, j, k, l]
            .into_iter()
            .circular_tuple_windows()
            .map(|(front, back)| {
                // Rather than search the edge list by floating-point endpoint
                // equality, use the canonical endpoint indices.
                let (edge_index, inverse) = EDGE_ENDPOINTS
                    .iter()
                    .enumerate()
                    .find_map(|(index, &[a, b])| {
                        if a == front && b == back {
                            Some((index, false))
                        } else if b == front && a == back {
                            Some((index, true))
                        } else {
                            None
                        }
                    })
                    .expect("face boundary must reference a canonical box edge");
                if inverse {
                    edges[edge_index].inverse()
                } else {
                    edges[edge_index].clone()
                }
            })
            .collect();
        Face::new_unchecked(vec![wire], bspline.into())
    };
    let shell: Shell = FACE_CORNERS
        .into_iter()
        .take(FACE_COUNT)
        .map(plane)
        .collect();
    (shell, edges, vertices)
}

/// Four-face semi-cube with splayed bottom vertices used for the open
/// two-edge fillet tests. This intentionally differs from the unit box.
fn build_stretched_open_box() -> (Shell, [Edge; 12], Vec<Vertex>) {
    let points = [
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
        Point3::new(0.0, -0.1, 0.0),
        Point3::new(1.1, -0.1, 0.0),
        Point3::new(1.1, 1.1, 0.0),
        Point3::new(0.0, 1.1, 0.0),
    ];
    build_box_from_points::<4>(points)
}

/// Add the two initial fillet side faces, then return the open top boundary
/// for the profile-specific second fillet. Keep original topology and edge
/// orientation, including which end edge is removed from the top wire.
fn prepare_stretched_open_wire(
    shell: &mut Shell,
    edges: &[Edge; 12],
    options: &FilletOptions,
) -> Wire {
    // Sequentially update both neighboring carrier faces and the shared
    // top face, preserving the original two-edge construction order.
    for (face0_index, face1_index, edge_index) in [(1, 2, 5), (2, 3, 6)] {
        let (face0, face1, face2, _, side1) = fillet_with_side(
            &shell[face0_index],
            &shell[face1_index],
            edges[edge_index].id(),
            None,
            Some(&shell[0]),
            options,
        )
        .unwrap();
        (shell[face0_index], shell[face1_index], shell[0]) = (face0, face1, side1.unwrap());
        shell.push(face2);
    }

    let mut boundary = shell[0].boundaries().pop().unwrap();
    boundary.pop_back();
    boundary
}

/// Exercise a five-face box's four-edge top loop with the chosen fillet
/// options. The profile-specific tests still supply their distinct geometry.
fn fillet_closed_top_loop(options: &FilletOptions) -> Shell {
    let (mut shell, edges, _) = build_5face_box();
    let face_count = shell.len();
    let top: Wire = edges[..4].iter().cloned().collect();
    assert!(top.is_closed());
    fillet_along_wire(&mut shell, &top, options).expect("closed top loop fillet");
    assert_eq!(shell.len(), face_count + 4);
    let _poly = shell.robust_triangulation(0.001).to_polygon();
    shell
}

/// Four-face test box: top/front/right/back, open on the bottom and left.
fn build_box_shell() -> (Shell, [Edge; 12], Vec<Vertex>) { build_box_with_faces::<4>() }

/// Five-face test box: top and four sides, open on the bottom.
fn build_5face_box() -> (Shell, [Edge; 12], Vec<Vertex>) { build_box_with_faces::<5>() }

/// Closed six-face test box.
fn build_6face_box() -> (Shell, [Edge; 12], Vec<Vertex>) { build_box_with_faces::<6>() }

#[test]
fn stretched_open_box_fixture_keeps_splayed_geometry_and_shared_edges() {
    let (mut shell, edges, vertices) = build_stretched_open_box();
    assert_eq!(shell.len(), 4);
    assert_eq!(vertices[4].point(), Point3::new(0.0, -0.1, 0.0));
    assert_eq!(vertices[5].point(), Point3::new(1.1, -0.1, 0.0));
    assert_eq!(vertices[6].point(), Point3::new(1.1, 1.1, 0.0));
    assert_eq!(vertices[7].point(), Point3::new(0.0, 1.1, 0.0));
    for (edge_index, first_face, second_face) in [(5, 1, 2), (6, 2, 3)] {
        assert!(
            shell[first_face]
                .edge_iter()
                .any(|e| e.id() == edges[edge_index].id())
        );
        assert!(
            shell[second_face]
                .edge_iter()
                .any(|e| e.id() == edges[edge_index].id())
        );
    }
    let opts = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        ..Default::default()
    };
    let boundary = prepare_stretched_open_wire(&mut shell, &edges, &opts);
    assert_eq!(shell.len(), 6);
    assert_eq!(boundary.front_vertex().unwrap(), &vertices[0]);
    assert!(!boundary.is_closed());
}
