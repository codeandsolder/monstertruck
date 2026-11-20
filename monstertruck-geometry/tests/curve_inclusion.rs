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

/// The original weighted rational plane fixtures, shared across curve
/// degrees to verify inclusion independently of surface weights.
fn rational_plane2() -> NurbsSurface<Vector3> {
    NurbsSurface::new(BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![
            vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.0, 2.0, 2.0)],
            vec![Vector3::new(2.0, 0.0, 2.0), Vector3::new(1.0, 1.0, 1.0)],
        ],
    ))
}

fn rational_plane3() -> NurbsSurface<Vector4> {
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
    NurbsSurface::new(BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        control_points,
    ))
}

#[test]
fn rational_surface_includes_weighted_curve_in_two_dimensions() {
    let surface = rational_plane2();
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
    let surface = rational_plane3();
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

/// A degree-zero spline has no derivatives but can jump between constant
/// pieces. Every span must be checked, not just the first control point.
#[test]
fn degree_zero_curve_inclusion_checks_all_spans_on_polynomial_surfaces() {
    let knots = KnotVector::from(vec![0.0, 1.0, 2.0, 3.0]);
    let inside2 = BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Point2::new(0.2, 0.2),
            Point2::new(0.6, 0.7),
            Point2::new(0.8, 0.4),
        ],
    );
    let outside2 = BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Point2::new(0.2, 0.2),
            Point2::new(1.4, 0.7),
            Point2::new(0.8, 0.4),
        ],
    );
    assert!(plane2().include(&inside2));
    assert!(!plane2().include(&outside2));

    let inside3 = BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Point3::new(0.2, 0.2, 0.0),
            Point3::new(0.6, 0.7, 0.0),
            Point3::new(0.8, 0.4, 0.0),
        ],
    );
    let outside3 = BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Point3::new(0.2, 0.2, 0.0),
            Point3::new(0.6, 0.7, 0.1),
            Point3::new(0.8, 0.4, 0.0),
        ],
    );
    assert!(plane3().include(&inside3));
    assert!(!plane3().include(&outside3));

    let inside_nurbs = NurbsCurve::new(BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Vector4::new(0.2, 0.2, 0.0, 1.0),
            Vector4::new(1.2, 1.4, 0.0, 2.0),
            Vector4::new(0.8, 0.4, 0.0, 1.0),
        ],
    ));
    let outside_nurbs = NurbsCurve::new(BsplineCurve::new_unchecked(
        knots,
        vec![
            Vector4::new(0.2, 0.2, 0.0, 1.0),
            Vector4::new(1.2, 1.4, 0.2, 2.0),
            Vector4::new(0.8, 0.4, 0.0, 1.0),
        ],
    ));
    assert!(plane3().include(&inside_nurbs));
    assert!(!plane3().include(&outside_nurbs));
}

/// The same shared inclusion sampler is used by rational 2D/3D surfaces.
#[test]
fn degree_zero_curve_inclusion_checks_all_spans_on_rational_surfaces() {
    let knots = KnotVector::from(vec![0.0, 1.0, 2.0, 3.0]);
    let surface2 = rational_plane2();
    let outside2 = NurbsCurve::new(BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Vector3::new(0.2, 0.2, 1.0),
            Vector3::new(1.4, 0.7, 1.0),
            Vector3::new(0.8, 0.4, 1.0),
        ],
    ));
    assert!(!surface2.include(&outside2));

    let surface3 = rational_plane3();
    let outside3 = BsplineCurve::new_unchecked(
        knots.clone(),
        vec![
            Point3::new(0.2, 0.2, 0.0),
            Point3::new(0.6, 0.7, 0.1),
            Point3::new(0.8, 0.4, 0.0),
        ],
    );
    let outside_nurbs = NurbsCurve::new(BsplineCurve::new_unchecked(
        knots,
        vec![
            Vector4::new(0.2, 0.2, 0.0, 1.0),
            Vector4::new(1.2, 1.4, 0.2, 2.0),
            Vector4::new(0.8, 0.4, 0.0, 1.0),
        ],
    ));
    assert!(!surface3.include(&outside3));
    assert!(!surface3.include(&outside_nurbs));
}
