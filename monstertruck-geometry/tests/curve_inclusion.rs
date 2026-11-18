//! Regression coverage for the six public B-spline/NURBS curve-inclusion
//! combinations. The shared sampling loop must preserve distinct polynomial
//! and rational surface projection strategies.

use monstertruck_geometry::prelude::*;

fn plane2() -> BsplineSurface<Point2> {
    BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![
            vec![Point2::new(0.0, 0.0), Point2::new(0.0, 1.0)],
            vec![Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)],
        ],
    )
}

fn plane3() -> BsplineSurface<Point3> {
    BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
    )
}

#[test]
fn polynomial_surface_includes_both_planar_curve_dimensions() {
    let curve2 = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point2::new(0.2, 0.4), Point2::new(0.8, 0.4)],
    );
    assert!(plane2().include(&curve2));
    let outside2 = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point2::new(0.2, 1.25), Point2::new(0.8, 1.25)],
    );
    assert!(!plane2().include(&outside2));

    let curve3 = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point3::new(0.2, 0.4, 0.0), Point3::new(0.8, 0.4, 0.0)],
    );
    assert!(plane3().include(&curve3));
    let off_surface3 = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point3::new(0.2, 0.4, 0.04), Point3::new(0.8, 0.4, 0.04)],
    );
    assert!(!plane3().include(&off_surface3));
}

#[test]
fn polynomial_surface_includes_weighted_rational_curve_only_on_plane() {
    let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![
            Vector4::new(0.2, 0.4, 0.0, 1.0),
            Vector4::new(1.6, 0.8, 0.0, 2.0),
        ],
    ));
    assert!(plane3().include(&curve));

    let off_surface = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![
            Vector4::new(0.2, 0.4, 0.02, 1.0),
            Vector4::new(1.6, 0.8, 0.04, 2.0),
        ],
    ));
    assert!(!plane3().include(&off_surface));
}

#[test]
fn rational_surface_includes_weighted_curve_in_two_dimensions() {
    let surface = NurbsSurface::new(BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![
            vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.0, 2.0, 2.0)],
            vec![Vector3::new(2.0, 0.0, 2.0), Vector3::new(1.0, 1.0, 1.0)],
        ],
    ));
    let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Vector3::new(0.2, 0.4, 1.0), Vector3::new(1.6, 0.8, 2.0)],
    ));
    assert!(surface.include(&curve));
    let outside = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Vector3::new(0.2, 1.4, 1.0), Vector3::new(1.6, 2.8, 2.0)],
    ));
    assert!(!surface.include(&outside));
}

#[test]
fn rational_surface_includes_polynomial_and_weighted_curves_in_three_dimensions() {
    let coords = [0.0_f64, 1.0];
    let control_points = coords
        .into_iter()
        .enumerate()
        .map(|(u_idx, u)| {
            coords
                .into_iter()
                .enumerate()
                .map(|(v_idx, v)| {
                    let weight = if u_idx == v_idx { 1.0 } else { 2.0 };
                    Vector4::new(u * weight, v * weight, 0.0, weight)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let surface = NurbsSurface::new(BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        control_points,
    ));
    let bspline = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point3::new(0.2, 0.4, 0.0), Point3::new(0.8, 0.4, 0.0)],
    );
    assert!(surface.include(&bspline));
    let outside_bspline = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point3::new(0.2, 0.4, 0.04), Point3::new(0.8, 0.4, 0.04)],
    );
    assert!(!surface.include(&outside_bspline));

    let nurbs = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![
            Vector4::new(0.2, 0.4, 0.0, 1.0),
            Vector4::new(1.6, 0.8, 0.0, 2.0),
        ],
    ));
    assert!(surface.include(&nurbs));
    let outside_nurbs = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![
            Vector4::new(0.2, 0.4, 0.04, 1.0),
            Vector4::new(1.6, 0.8, 0.08, 2.0),
        ],
    ));
    assert!(!surface.include(&outside_nurbs));
}
