//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;
use monstertruck_geometry::prelude::TryIntoHomogeneousBsplineSurface;

/// The bound around a ball is the ball's own box, and it is looser than the
/// ball by `6 / pi` -- 91.0%. Stated against the closed form so the number
/// is checked, not asserted.
#[test]
fn an_analytic_sphere_is_bounded_by_its_own_box_and_by_6_over_pi() {
    let center = Point3::new(3.0, -4.0, 0.5);
    let radius = 12.5;
    let surface = Surface::SphericalSurface(Processor::new(Sphere::new(center, radius)));
    let bound = certified_surface_bounding_box(&surface, None)
        .expect("an analytic sphere bounds itself without a boundary hull");
    let diagonal = bound.diagonal();
    for (axis, extent) in [("x", diagonal.x), ("y", diagonal.y), ("z", diagonal.z)] {
        assert!(
            (extent - 2.0 * radius).abs() <= 1.0e-12 * radius,
            "a ball of radius {radius} must box to {} along {axis}; got {extent}",
            2.0 * radius,
        );
    }
    let box_volume = diagonal.x * diagonal.y * diagonal.z;
    let ball_volume = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    let looseness = box_volume / ball_volume;
    assert!(
        (looseness - ANALYTIC_SPHERE_LOOSENESS).abs() <= 1.0e-12,
        "the certified box around a ball must exceed it by exactly 6/pi = \
         {ANALYTIC_SPHERE_LOOSENESS}; measured {looseness}",
    );
    // Every point of the sphere is inside, including the six poles that
    // carry no vertex -- the failure mode the vertex box has.
    for offset in [
        Vector3::unit_x(),
        -Vector3::unit_x(),
        Vector3::unit_y(),
        -Vector3::unit_y(),
        Vector3::unit_z(),
        -Vector3::unit_z(),
    ] {
        assert!(
            bound.contains(center + offset * radius),
            "the pole {:?} of the ball must be inside its certified box",
            center + offset * radius,
        );
    }
}

/// The SAME ball routed the other way -- as the rational NURBS net the
/// repo's own `TryIntoHomogeneousBsplineSurface` builds for a `Sphere`.
///
/// A certified control-hull box is allowed to be conservative. Its correctness
/// contract is enclosure, not reproducing the current control net's exact box.
#[test]
fn rational_sphere_control_hull_certifies_the_entire_sphere() {
    let radius = 2.0;
    let center = Point3::new(-1.0, 0.5, 4.0);
    let sphere = Sphere::new(center, radius);
    let net = sphere
        .try_into_homogeneous_bspline_surface()
        .expect("the repo builds a rational net for a sphere");
    let surface = Surface::NurbsSurface(NurbsSurface::new(net));
    let bound = certified_surface_bounding_box(&surface, None)
        .expect("a positively-weighted rational patch bounds itself");

    // Check the actual geometric surface, not a captured property of today's
    // NURBS control net. The poles are included explicitly by the latitude grid.
    for latitude_step in 0..=32 {
        let latitude =
            std::f64::consts::PI * f64::from(latitude_step) / 32.0 - std::f64::consts::FRAC_PI_2;
        for longitude_step in 0..64 {
            let longitude = std::f64::consts::TAU * f64::from(longitude_step) / 64.0;
            let point = center
                + Vector3::new(
                    radius * latitude.cos() * longitude.cos(),
                    radius * latitude.cos() * longitude.sin(),
                    radius * latitude.sin(),
                );
            assert!(
                bound.contains(point),
                "sphere point {point:?} escapes the rational net's certified control hull",
            );
        }
    }
}

/// Composite curve bounds must come from structural carrier evidence, not the
/// fitted leader. Two spheres intersect in the unit circle at z = 0; their
/// certified carrier boxes alone therefore enclose the represented curve.
#[test]
fn intersection_curve_exposes_its_structural_certified_bound() {
    let radius = f64::sqrt(2.0);
    let surface0 = Surface::SphericalSurface(Processor::new(Sphere::new(
        Point3::new(0.0, 0.0, 1.0),
        radius,
    )));
    let surface1 = Surface::SphericalSurface(Processor::new(Sphere::new(
        Point3::new(0.0, 0.0, -1.0),
        radius,
    )));
    let leader = Curve::NurbsCurve(
        TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, std::f64::consts::TAU))
            .to_same_geometry(),
    );
    let curve = Curve::IntersectionCurve(SurfaceCurve::with_boundaries(
        Box::new(surface0),
        Box::new(surface1),
        Box::new(leader),
        None,
        None,
    ));

    let direct = certified_curve_bounding_box(&curve)
        .expect("the two self-bounded carrier surfaces certify the intersection curve");
    let via_trait = curve
        .certified_spatial_bounds()
        .expect("BoundedCurve must expose the same certified structural enclosure");
    assert_eq!(direct.min(), via_trait.min());
    assert_eq!(direct.max(), via_trait.max());

    // The exact carrier intersection is the unit circle in z = 0. Check the
    // structural enclosure against that geometry directly instead of depending
    // on the iterative intersection evaluator, which has independent singular
    // convergence behavior at some parameters.
    for step in 0..64 {
        let angle = std::f64::consts::TAU * f64::from(step) / 64.0;
        let point = Point3::new(angle.cos(), angle.sin(), 0.0);
        assert!(
            direct.contains(point),
            "exact sphere-sphere intersection point {point:?} escapes its certified carrier box",
        );
    }
}

/// A non-positive weight voids the rational convex-hull property, so the
/// bound refuses rather than inventing one.
#[test]
fn a_non_positive_weight_refuses() {
    let net = vec![
        vec![
            Vector4::new(0.0, 0.0, 0.0, 1.0),
            Vector4::new(0.0, 1.0, 0.0, 1.0),
        ],
        vec![
            Vector4::new(1.0, 0.0, 0.0, 1.0),
            // The weight that voids the property.
            Vector4::new(1.0, 1.0, 0.0, 0.0),
        ],
    ];
    let knots = KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]);
    let surface = Surface::NurbsSurface(NurbsSurface::new(BsplineSurface::new_unchecked(
        (knots.clone(), knots),
        net,
    )));
    assert!(
        certified_surface_bounding_box(&surface, None).is_none(),
        "a zero weight must refuse: the convex-hull property does not hold",
    );
}

/// A cylinder's certified box is the tube its own boundary circles cut out
/// -- not the infinite surface, and not the vertex hull.
#[test]
fn a_cylinder_face_is_bounded_by_its_boundary_wires() {
    let radius = 3.0;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let axis = Vector3::unit_z();
    let profile = Line(Point3::new(radius, 0.0, 0.0), Point3::new(radius, 0.0, 1.0));
    let surface = Surface::RevolutionSurface(Processor::new(RevolutionSurface::by_revolution(
        Curve::Line(profile),
        origin,
        axis,
    )));
    // Two boundary circles at z = 0 and z = 7, given as their (sound)
    // enclosing squares -- a rational circle's control hull is exactly that.
    let hull: Vec<Point3> = [0.0, 7.0]
        .into_iter()
        .flat_map(|z| {
            [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
                .into_iter()
                .map(move |(x, y)| Point3::new(x * radius, y * radius, z))
        })
        .collect();
    let bound = certified_surface_bounding_box(&surface, Some(&hull))
        .expect("a straight-profile revolution bounds from its wires");
    let diagonal = bound.diagonal();
    // The hull's corners sit at radius*sqrt(2) from the axis, so the tube is
    // bounded at that radius -- sound, and looser than the true 3.0.
    let expected = 2.0 * radius * std::f64::consts::SQRT_2;
    assert!(
        (diagonal.x - expected).abs() <= 1.0e-12 * expected
            && (diagonal.y - expected).abs() <= 1.0e-12 * expected,
        "the tube's radial extent must be the hull's {expected}; got {} x {}",
        diagonal.x,
        diagonal.y,
    );
    assert!(
        (diagonal.z - 7.0).abs() <= 1.0e-12,
        "the tube's axial extent is exactly its wires' 7.0; got {}",
        diagonal.z,
    );
    // Every point of the cylinder wall is inside.
    for step in 0..16 {
        let angle = std::f64::consts::TAU * f64::from(step) / 16.0;
        for z in [0.0, 3.5, 7.0] {
            let point = Point3::new(radius * angle.cos(), radius * angle.sin(), z);
            assert!(
                bound.contains(point),
                "the wall point {point:?} must be inside the certified box",
            );
        }
    }
}

/// Without a boundary hull a cylinder cannot be bounded at all, and says so.
#[test]
fn a_cylinder_without_wires_refuses() {
    let surface = Surface::RevolutionSurface(Processor::new(RevolutionSurface::by_revolution(
        Curve::Line(Line(Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0))),
        Point3::origin(),
        Vector3::unit_z(),
    )));
    assert!(
        certified_surface_bounding_box(&surface, None).is_none(),
        "an untrimmed cylinder is unbounded; a bound must not be invented",
    );
}

/// A ring torus is bounded by `R + r` across and `r` along its axis.
#[test]
fn a_ring_torus_is_bounded_by_its_own_radii() {
    let (large, small) = (5.0, 1.5);
    let center = Point3::new(1.0, 2.0, 3.0);
    let surface =
        Surface::ToroidalSurface(Processor::new(Torus::new_unchecked(center, large, small)));
    let bound = certified_surface_bounding_box(&surface, None).expect("a torus bounds itself");
    let diagonal = bound.diagonal();
    assert!(
        (diagonal.x - 2.0 * (large + small)).abs() <= 1.0e-12
            && (diagonal.y - 2.0 * (large + small)).abs() <= 1.0e-12,
        "a ring torus spans 2(R+r) across its axis; got {} x {}",
        diagonal.x,
        diagonal.y,
    );
    assert!(
        (diagonal.z - 2.0 * small).abs() <= 1.0e-12,
        "a ring torus spans 2r along its axis; got {}",
        diagonal.z,
    );
}

/// The whole point, in one assertion: on the `#25387` geometry the vertex
/// box is SMALLER than the solid's own volume and the certified box is not.
///
/// The numbers are the measured ones (spec 013): `R = 12.5`, `h = 9`,
/// `r = 6`, closed-form volume `5273.16`, vertex box `18 x 17.349 x 12`.
/// Only the sphere face is needed to make the point, so this stays a pure
/// geometry test with no corpus.
#[test]
fn the_certified_box_contains_a_volume_the_vertex_box_does_not() {
    let (big_radius, half_x, bore_radius) = (12.5f64, 9.0f64, 6.0f64);
    let trim = (big_radius * big_radius - half_x * half_x).sqrt();
    let closed_form = 2.0
        * std::f64::consts::PI
        * (big_radius * big_radius * half_x
            - half_x.powi(3) / 3.0
            - bore_radius * bore_radius * half_x);
    let vertex_box = (2.0 * half_x) * (2.0 * trim) * (2.0 * bore_radius);
    assert!(
        closed_form > vertex_box,
        "the witness requires the correct volume {closed_form} to EXCEED the \
         vertex box {vertex_box}",
    );
    let surface =
        Surface::SphericalSurface(Processor::new(Sphere::new(Point3::origin(), big_radius)));
    let bound = certified_surface_bounding_box(&surface, None).expect("sphere bounds itself");
    let diagonal = bound.diagonal();
    let certified = diagonal.x * diagonal.y * diagonal.z;
    assert!(
        certified > closed_form,
        "the certified box {certified} must contain the solid's volume \
         {closed_form}; that is the whole defect",
    );
}
