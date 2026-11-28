//! A trim that is closed in space but open in parameter space -- a
//! cylinder's circle, running `u` from 0 to 2π -- matches a point at its
//! start with both of its ends. Projecting that point without a hint must
//! still land on the start: a boundary is walked from there, and landing on
//! the far end sends every later projection off the curve and folds the
//! face.

use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI, TAU};

use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;

/// A rational quadratic cylinder of radius 1 and height 1 around `z`, `u`
/// over `[0, 2π]` in four quarter spans, `v` over `[0, 1]`.
fn cylinder() -> NurbsSurface<Vector4> {
    let ring = [
        (1.0, 0.0, 1.0),
        (1.0, 1.0, FRAC_1_SQRT_2),
        (0.0, 1.0, 1.0),
        (-1.0, 1.0, FRAC_1_SQRT_2),
        (-1.0, 0.0, 1.0),
        (-1.0, -1.0, FRAC_1_SQRT_2),
        (0.0, -1.0, 1.0),
        (1.0, -1.0, FRAC_1_SQRT_2),
        (1.0, 0.0, 1.0),
    ];
    let control_points = ring
        .iter()
        .map(|&(x, y, w)| {
            [0.0, 1.0]
                .iter()
                .map(|&z| Vector4::new(x * w, y * w, z * w, w))
                .collect()
        })
        .collect();
    let u_knots = KnotVector::from(vec![
        0.0,
        0.0,
        0.0,
        FRAC_PI_2,
        FRAC_PI_2,
        PI,
        PI,
        3.0 * FRAC_PI_2,
        3.0 * FRAC_PI_2,
        TAU,
        TAU,
        TAU,
    ]);
    let v_knots = KnotVector::bezier_knot(1);
    NurbsSurface::new(BsplineSurface::new_unchecked(
        (u_knots, v_knots),
        control_points,
    ))
}

/// Projects both evaluations of `trim`'s meeting point -- a neighbouring
/// face's sample of it can be either, a rounding apart -- and requires each
/// to land at the trim's start.
fn assert_projects_to_start(trim: &ParameterCurve<BsplineCurve<Point2>, NurbsSurface<Vector4>>) {
    let (start, end) = trim.range_tuple();
    assert!(
        trim.evaluate(start).near(&trim.evaluate(end)),
        "the fixture must be closed in space"
    );
    for at in [start, end] {
        let (parameter, uv) = trim
            .project_boundary_point(trim.evaluate(at), None, KERNEL_EPSILON)
            .expect("the start of the trim projects onto it");
        assert!(
            (parameter - start).abs() < (parameter - end).abs(),
            "the point at {at} projected to {parameter} (uv {uv:?}), \
             the far end of [{start}, {end}]"
        );
    }
}

/// A circle's trim, from `(from, v)` to `(to, v)`.
fn circle_trim(
    from: f64,
    to: f64,
    v: f64,
) -> ParameterCurve<BsplineCurve<Point2>, NurbsSurface<Vector4>> {
    let line = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point2::new(from, v), Point2::new(to, v)],
    );
    ParameterCurve::new(line, cylinder())
}

#[test]
fn a_trim_running_up_the_seam_projects_its_start_to_its_start() {
    assert_projects_to_start(&circle_trim(0.0, TAU, 0.0));
}

/// The case a real part hit: the top circle of a hole, taken in its trim
/// loop's direction, runs `u` down from 2π.
#[test]
fn a_trim_running_down_the_seam_projects_its_start_to_its_start() {
    assert_projects_to_start(&circle_trim(TAU, 0.0, 1.0));
}

#[test]
fn refuses_extrapolated_trim_projection() {
    let trim = ParameterCurve::new(
        Line(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
        Plane::xy(),
    );
    assert!(
        trim.project_boundary_point(Point3::new(2.0, 0.0, 0.0), None, 1e-6)
            .is_none()
    );
    assert!(
        trim.project_boundary_point(Point3::new(0.5, 0.2, 0.0), None, 1e-6)
            .is_none()
    );
    let (t, uv) = trim
        .project_boundary_point(Point3::new(0.5, 0.0, 0.0), Some(1e12), 1e-6)
        .expect("in-range projection");
    assert!((0.0..=1.0).contains(&t));
    assert!(uv.near(&Point2::new(0.5, 0.0)));
}

#[test]
fn seam_crossing_rational_face_consumes_every_shared_boundary_sample() {
    use monstertruck_modeling::{Curve, Curve2D, Surface};
    use monstertruck_topology::{
        compress::{
            CompressedEdge, CompressedEdgeUse, CompressedTrimmedFace, CompressedTrimmedShell,
        },
        shell::ShellCondition,
    };
    type Trim = ParameterCurve<Curve2D, Box<Surface>>;
    // The parameter curve has no supplied UV trim. Its rational spatial
    // arc crosses the endpoint of the finite clamped circle chart.
    let theta = 0.35_f64;
    let c = theta.cos();
    let s = theta.sin();
    let vertices = vec![
        Point3::new(c, -s, 0.0),
        Point3::new(c, s, 0.0),
        Point3::new(c, s, 1.0),
        Point3::new(c, -s, 1.0),
    ];
    let arc = |z: f64| {
        Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(c, -s, z, 1.0),
                Vector4::new(1.0, 0.0, z * c, c),
                Vector4::new(c, s, z, 1.0),
            ],
        )))
    };
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: arc(0.0),
        },
        CompressedEdge {
            vertices: (1, 2),
            curve: Curve::Line(Line(vertices[1], vertices[2])),
        },
        CompressedEdge {
            vertices: (3, 2),
            curve: arc(1.0),
        },
        CompressedEdge {
            vertices: (3, 0),
            curve: Curve::Line(Line(vertices[3], vertices[0])),
        },
    ];
    let shell = CompressedTrimmedShell::<Point3, Curve, Surface, Trim> {
        vertices,
        edges,
        faces: vec![CompressedTrimmedFace {
            boundaries: vec![
                [(0, true), (1, true), (2, false), (3, true)]
                    .into_iter()
                    .map(|(index, orientation)| CompressedEdgeUse {
                        index,
                        orientation,
                        trim_curve: None,
                    })
                    .collect(),
            ],
            orientation: true,
            surface: Surface::NurbsSurface(cylinder()),
        }],
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };
    let tolerance = 0.001;
    let meshed = shell.triangulation(tolerance);
    let mesh = stitch_compressed_shell(&meshed, KERNEL_EPSILON)
        .expect("consume exact native sample chains");
    assert_eq!(mesh.shell_condition(), ShellCondition::Oriented);
    let area = mesh
        .face_iter()
        .map(|polygon| {
            let origin = mesh.positions()[polygon[0].pos];
            polygon[1..]
                .windows(2)
                .map(|pair| {
                    let a = mesh.positions()[pair[0].pos] - origin;
                    let b = mesh.positions()[pair[1].pos] - origin;
                    a.cross(&b).norm() * 0.5
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    assert!(
        2.0_f64.mul_add(-theta, area).abs() < 0.005,
        "area {area} instead of {}",
        2.0 * theta
    );
    assert!(
        mesh.positions()
            .iter()
            .all(|p| (p.x.hypot(p.y) - 1.0).abs() < KERNEL_EPSILON)
    );
}
