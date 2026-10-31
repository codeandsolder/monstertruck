//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;

fn plane_surface() -> Surface {
    Surface::Plane(Plane::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ))
}

fn exact_quadratic_leader() -> NurbsCurve<Vector4> {
    NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(2),
        vec![
            Vector4::new(0.0, 0.0, 0.0, 1.0),
            Vector4::new(0.5, 1.0, 0.0, 1.0),
            Vector4::new(1.0, 0.0, 0.0, 1.0),
        ],
    ))
}

fn line_pcurve_on_plane() -> Curve {
    Curve::ParameterCurve(ParameterCurve::new(
        Curve2D::Line(Line(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0))),
        Box::new(plane_surface()),
    ))
}

fn intersection_curve_with_leader(leader: Curve) -> Curve {
    Curve::Intersection(SurfaceCurve::with_boundaries(
        Box::new(plane_surface()),
        Box::new(plane_surface()),
        Box::new(leader),
        None,
        None,
    ))
}

/// Boolean seam edges are `IntersectionCurve`s with exact leaders. The
/// fillet conversion must return THAT curve, not a resampled polyline.
#[test]
fn intersection_curve_with_exact_leader_converts_exactly() {
    let leader = exact_quadratic_leader();
    let curve = intersection_curve_with_leader(Curve::NurbsCurve(leader.clone()));
    let converted = NurbsCurve::<Vector4>::try_from(curve).expect("exact leader must convert");
    assert_eq!(converted.degree(), leader.degree());
    assert_eq!(converted.knot_vector(), leader.knot_vector());
    assert_eq!(converted.control_points(), leader.control_points());
}

/// A parameter curve has no exact NURBS representation: refuse, never
/// silently sample.
#[test]
fn parameter_curve_refuses_instead_of_sampling() {
    assert!(NurbsCurve::<Vector4>::try_from(line_pcurve_on_plane()).is_err());
}

/// A seam whose leader chain bottoms out non-exact must refuse as well.
#[test]
fn intersection_curve_with_non_exact_leader_refuses() {
    let curve = intersection_curve_with_leader(line_pcurve_on_plane());
    assert!(NurbsCurve::<Vector4>::try_from(curve).is_err());
}

/// The generic modeling bridge must preserve round-fillet intersection curves
/// exactly enough for the rebuilt modeling shell to remain geometrically consistent.
#[test]
fn generic_round_fillet_preserves_modeling_consistency() {
    use crate::builder;
    use monstertruck_fillet::{FilletOptions, fillet_edges_generic};

    let vertex = builder::vertex(Point3::origin());
    let edge: monstertruck_topology::Edge<Point3, Curve> =
        builder::extrude(&vertex, Vector3::new(4.0, 0.0, 0.0));
    let face: monstertruck_topology::Face<Point3, Curve, Surface> =
        builder::extrude(&edge, Vector3::new(0.0, 3.0, 0.0));
    let solid: monstertruck_topology::Solid<Point3, Curve, Surface> =
        builder::extrude(&face, Vector3::new(0.0, 0.0, 2.0));
    let mut boundaries = solid.into_boundaries();
    let shell = boundaries.first_mut().expect("box must have one shell");
    let start = Point3::new(0.0, 0.0, 0.0);
    let end = Point3::new(4.0, 0.0, 0.0);
    let selected = shell
        .edge_iter()
        .find(|edge| {
            let front = edge.absolute_front().point();
            let back = edge.absolute_back().point();
            (front.near(&start) && back.near(&end)) || (front.near(&end) && back.near(&start))
        })
        .expect("box edge must exist")
        .clone();

    fillet_edges_generic(shell, &[selected], Some(&FilletOptions::constant(0.25)))
        .expect("round fillet should succeed");
    assert!(shell.is_geometric_consistent());
}
