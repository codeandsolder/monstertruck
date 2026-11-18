//! Baseline behavioral contracts for STEP conic trimming, both dimensions.
//! In particular, STEP parabolas must remain parabolas rather than being
//! accidentally reinterpreted as hyperbolas.
use super::super::{
    Axis2Placement, Axis2Placement2d, Axis2Placement3d, Circle, Ellipse, Hyperbola, Parabola,
};
use super::*;
use geom::{BoundedCurve, ParametricCurve, Tolerance};

fn placement(dim: usize) -> Axis2Placement {
    let location = CartesianPoint {
        label: String::new(),
        coordinates: if dim == 2 {
            vec![0.0, 0.0]
        } else {
            vec![0.0, 0.0, 0.0]
        },
    };
    if dim == 2 {
        Axis2Placement::Axis2Placement2d(Axis2Placement2d {
            label: String::new(),
            location,
            ref_direction: None,
        })
    } else {
        Axis2Placement::Axis2Placement3d(Axis2Placement3d {
            label: String::new(),
            location,
            axis: None,
            ref_direction: None,
        })
    }
}

fn conic(dim: usize, kind: &str, scale: f64) -> Conic {
    let position = placement(dim);
    let label = String::new();
    match kind {
        "circle" => Conic::Circle(Circle {
            label,
            position,
            radius: 2.0 * scale,
        }),
        "ellipse" => Conic::Ellipse(Ellipse {
            label,
            position,
            semi_axis_1: 2.0 * scale,
            semi_axis_2: 3.0 * scale,
        }),
        "hyperbola" => Conic::Hyperbola(Hyperbola {
            label,
            position,
            semi_axis: 2.0 * scale,
            semi_imag_axis: 3.0 * scale,
        }),
        "parabola" => Conic::Parabola(Parabola {
            label,
            position,
            focal_dist: 2.0 * scale,
        }),
        _ => unreachable!(),
    }
}

fn point(kind: &str, t: f64) -> Point2 {
    match kind {
        "circle" => Point2::new(2.0 * t.cos(), 2.0 * t.sin()),
        "ellipse" => Point2::new(2.0 * t.cos(), 3.0 * t.sin()),
        "hyperbola" => Point2::new(2.0 * t.cosh(), 3.0 * t.sinh()),
        "parabola" => Point2::new(2.0 * t * t, 4.0 * t),
        _ => unreachable!(),
    }
}

macro_rules! check_conic_endpoints {
    ($parsed:expr, $a:expr, $b:expr, $kind:expr, $dim:literal, $curve_ty:ident, $conic_ty:ident) => {{
        let parsed = $parsed.unwrap_or_else(|err| panic!("{} {}: {err}", $dim, $kind));
        let (start, end) = BoundedCurve::range_tuple(&parsed);
        assert!(
            ParametricCurve::evaluate(&parsed, start).near(&$a),
            "{} {} start",
            $dim,
            $kind
        );
        assert!(
            ParametricCurve::evaluate(&parsed, end).near(&$b),
            "{} {} end",
            $dim,
            $kind
        );
        assert!(
            matches!(
                ($kind, parsed),
                (
                    "circle" | "ellipse",
                    $curve_ty::Conic($conic_ty::Ellipse(_))
                ) | ("hyperbola", $curve_ty::Conic($conic_ty::Hyperbola(_)))
                    | ("parabola", $curve_ty::Conic($conic_ty::Parabola(_)))
            ),
            "unexpected conic variant {} {}",
            $dim,
            $kind,
        );
    }};
}

#[test]
fn all_conic_variants_trim_2d_and_3d_without_moving_endpoints() {
    for kind in ["circle", "ellipse", "hyperbola", "parabola"] {
        let a = point(kind, 0.2);
        let b = point(kind, 0.8);
        check_conic_endpoints!(
            EdgeCurve::parse_conic_2d(&conic(2, kind, 1.0), a, b),
            a,
            b,
            kind,
            "2D",
            Curve2D,
            Conic2D
        );
        let a3 = Point3::new(a.x, a.y, 0.0);
        let b3 = Point3::new(b.x, b.y, 0.0);
        check_conic_endpoints!(
            EdgeCurve::parse_conic_3d(&conic(3, kind, 1.0), a3, b3),
            a3,
            b3,
            kind,
            "3D",
            Curve3D,
            Conic3D
        );
    }
}

#[test]
fn singular_conic_transform_errors_remain_dimension_specific() {
    for kind in ["circle", "ellipse", "hyperbola", "parabola"] {
        let p2 = Point2::new(1.0, 0.0);
        let p3 = Point3::new(1.0, 0.0, 0.0);
        let err2 = EdgeCurve::parse_conic_2d(&conic(2, kind, 0.0), p2, p2)
            .unwrap_err()
            .to_string();
        let err3 = EdgeCurve::parse_conic_3d(&conic(3, kind, 0.0), p3, p3)
            .unwrap_err()
            .to_string();
        let expected2 = match kind {
            "circle" | "ellipse" => "Failed to convert Circle",
            "hyperbola" => "Failed to convert Hyperbola",
            _ => "Failed to convert Parabola",
        };
        let expected3 = match kind {
            "circle" | "ellipse" | "hyperbola" => "Failed to convert Circle",
            _ => "Failed to convert Parabola",
        };
        assert_eq!(err2, expected2, "2D {kind}");
        assert_eq!(err3, expected3, "3D {kind}");
    }
}

#[test]
fn hyperbolas_do_not_wrap_nonperiodic_parameters() {
    let p = point("hyperbola", 0.6);
    let q = point("hyperbola", -0.2);
    let curve2 = EdgeCurve::parse_conic_2d(&conic(2, "hyperbola", 1.0), p, q).unwrap();
    let curve3 = EdgeCurve::parse_conic_3d(
        &conic(3, "hyperbola", 1.0),
        Point3::new(p.x, p.y, 0.0),
        Point3::new(q.x, q.y, 0.0),
    )
    .unwrap();
    let (start2, end2) = BoundedCurve::range_tuple(&curve2);
    let (start3, end3) = BoundedCurve::range_tuple(&curve3);
    assert!(
        end2 < start2,
        "2D reversed hyperbola trim must be decreasing"
    );
    assert!(
        end3 < start3,
        "3D reversed hyperbola trim must be decreasing"
    );
    assert!((start2 - start3).abs() < 1.0e-8);
    assert!((end2 - end3).abs() < 1.0e-8);
    assert!(ParametricCurve::evaluate(&curve2, start2).near(&p));
    assert!(ParametricCurve::evaluate(&curve2, end2).near(&q));
    assert!(ParametricCurve::evaluate(&curve3, start3).near(&Point3::new(p.x, p.y, 0.0)));
    assert!(ParametricCurve::evaluate(&curve3, end3).near(&Point3::new(q.x, q.y, 0.0)));
}
