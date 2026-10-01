//! Downstream-application smoke scenario: build primitive solids, combine them
//! with the boolean operators, and confirm the results tessellate into
//! non-degenerate watertight meshes.
//!
//! This mirrors how an interactive CAD app (e.g. an OpenCADStudio-style scene
//! graph) drives the kernel: a handful of primitives fed through
//! `union`/`intersection`/`difference`, then meshed for display. It is a
//! coarse "does the pipeline stay closed end to end" check, not a numerical
//! accuracy test -- the `transversal::integrate` unit tests own the exact
//! values.
//!
//! Note: `monstertruck_solid::{and, or}` build their output through
//! `Solid::try_new`, which only returns `Ok` when every boundary shell is a
//! closed manifold. So a returned `Solid` is watertight by construction; these
//! tests additionally confirm it tessellates to a non-empty mesh whose volume
//! matches the analytic result. Because every primitive here is flat-faced, the
//! tessellated volume is exact (independent of mesh density), so an exact-value
//! assertion is both legitimate and a far stronger regression guard than a
//! "non-empty output" smoke check -- it would catch a boolean that silently
//! kept the wrong region.

use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::*;

const TOL: f64 = 0.05;
const VOLUME_EPS: f64 = 1.0e-3;

/// Axis-aligned unit cube with its minimum corner at `origin`.
fn unit_cube(origin: Point3) -> Solid {
    let v = builder::vertex(origin);
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    builder::extrude(&f, Vector3::unit_z())
}

/// Square-section vertical column centered on `(cx, cy)`, tall enough to pierce
/// a unit cube standing on `z = 0`.
fn square_column(cx: f64, cy: f64, half: f64) -> Solid {
    let v = builder::vertex(Point3::new(cx - half, cy - half, -0.5));
    let e = builder::extrude(&v, Vector3::unit_x() * (2.0 * half));
    let f = builder::extrude(&e, Vector3::unit_y() * (2.0 * half));
    builder::extrude(&f, Vector3::unit_z() * 2.0)
}

/// Asserts a boolean result is watertight and geometrically correct: it has at
/// least one boundary shell, tessellates to a non-empty mesh, and that mesh has
/// the expected (positively oriented) volume.
fn assert_solid(label: &str, solid: &Solid, expected_volume: f64) -> Result<()> {
    anyhow::ensure!(
        !solid.boundaries().is_empty(),
        "{label}: solid has no boundary shells"
    );
    let mesh = solid.triangulation(0.01).to_polygon();
    let triangles = mesh.faces().triangle_iter().count();
    anyhow::ensure!(triangles > 0, "{label}: tessellation produced no triangles");
    let volume = mesh.volume();
    anyhow::ensure!(
        (volume - expected_volume).abs() < VOLUME_EPS,
        "{label}: tessellated volume {volume:.6}, expected {expected_volume:.6}"
    );
    Ok(())
}

#[test]
fn union_of_overlapping_cubes() -> Result<()> {
    // Two unit cubes overlapping in a `0.5` cube: `2 - 0.5^3 = 1.875`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::or(&a, &b, TOL)?;
    assert_solid("union of overlapping cubes", &result, 1.875)
}

#[test]
fn intersection_of_overlapping_cubes() -> Result<()> {
    // The shared `0.5` cube: `0.5^3 = 0.125`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::and(&a, &b, TOL)?;
    assert_solid("intersection of overlapping cubes", &result, 0.125)
}

#[test]
fn difference_of_overlapping_cubes() -> Result<()> {
    // First cube minus the shared region: `1 - 0.125 = 0.875`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::difference(&a, &b, TOL)?;
    assert_solid("difference of overlapping cubes", &result, 0.875)
}

#[test]
fn cube_minus_column() -> Result<()> {
    // Unit cube with a `0.4 x 0.4` column punched through: `1 - 0.4^2 = 0.84`.
    let cube = unit_cube(Point3::origin());
    let column = square_column(0.5, 0.5, 0.2);
    let result = monstertruck_solid::difference(&cube, &column, TOL)?;
    assert_solid("cube minus square column", &result, 0.84)
}

#[test]
fn profile_generated_difference_matrix() -> Result<()> {
    fn profile_box(min: Point3, max: Point3) -> Result<Solid> {
        let v = builder::vertices([
            Point3::new(min.x, min.y, min.z),
            Point3::new(max.x, min.y, min.z),
            Point3::new(max.x, max.y, min.z),
            Point3::new(min.x, max.y, min.z),
        ]);
        let wire: Wire = vec![
            builder::line(&v[0], &v[1]),
            builder::line(&v[1], &v[2]),
            builder::line(&v[2], &v[3]),
            builder::line(&v[3], &v[0]),
        ]
        .into();
        profile::solid_from_planar_profile::<Curve, Surface>(
            vec![wire],
            Vector3::new(0.0, 0.0, max.z - min.z),
        )
        .map_err(Into::into)
    }

    let primitive_host: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 6.0, 2.0),
    ]));
    let primitive_cutter: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(3.0, 2.0, -1.0),
        Point3::new(7.0, 4.0, 3.0),
    ]));
    let profile_host = profile_box(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 6.0, 2.0),
    )?;
    let profile_cutter = profile_box(
        Point3::new(3.0, 2.0, -1.0),
        Point3::new(7.0, 4.0, 3.0),
    )?;
    let profile_cutter_at_origin = profile_box(
        Point3::new(3.0, 2.0, 0.0),
        Point3::new(7.0, 4.0, 4.0),
    )?;
    let transformed_profile_cutter = builder::transformed(
        &profile_cutter_at_origin,
        Matrix4::from_translation(Vector3::new(0.0, 0.0, -1.0)),
    );

    for (label, host, cutter) in [
        ("primitive/primitive", &primitive_host, &primitive_cutter),
        ("profile/primitive", &profile_host, &primitive_cutter),
        ("primitive/profile", &primitive_host, &profile_cutter),
        ("profile/profile", &profile_host, &profile_cutter),
        (
            "profile/transformed-profile",
            &profile_host,
            &transformed_profile_cutter,
        ),
    ] {
        match monstertruck_solid::difference(host, cutter, TOL) {
            Ok(result) => {
                eprintln!("{label}: ok");
                assert_solid(label, &result, 104.0)?;
            }
            Err(error) => eprintln!("{label}: {error:?}"),
        }
    }

    // Keep the downstream construction as the actual regression assertion.
    let result = monstertruck_solid::difference(
        &profile_host,
        &transformed_profile_cutter,
        TOL,
    )?;
    assert_solid(
        "profile-generated box minus transformed column",
        &result,
        104.0,
    )
}

#[test]
fn cube_minus_enclosed_cube_preserves_inward_cavity_orientation() -> Result<()> {
    let outer = unit_cube(Point3::origin());
    let inner: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.25, 0.25, 0.25),
        Point3::new(0.75, 0.75, 0.75),
    ]));
    let result = monstertruck_solid::difference(&outer, &inner, TOL)?;
    anyhow::ensure!(
        result.boundaries().len() == 2,
        "enclosed subtraction must produce outer and cavity boundary shells"
    );
    assert_solid("cube minus enclosed cube", &result, 0.875)
}

/// Regression for a radial slot cut from a solid of revolution.
///
/// The cutter reaches the revolution axis and exits through the top end. This
/// is the minimal form of the failure reported by downstream STEP recovery:
/// the classic boolean pipeline used to classify every divided face out of the
/// AND result and return `EmptyOutputShell`.
#[test]
fn revolved_cylinder_minus_axis_touching_radial_slot() -> Result<()> {
    let vertices = builder::vertices([
        Point3::new(0.0, 0.0, -2.0),
        Point3::new(2.0, 0.0, -2.0),
        Point3::new(2.0, 0.0, 2.0),
        Point3::new(0.0, 0.0, 2.0),
    ]);
    // Use the axis-aware wire revolution path. Generic face revolution would
    // sweep the closing on-axis edge into zero-length topology, which is not a
    // valid input for testing the Boolean itself.
    let profile: Wire = vec![
        builder::line(&vertices[0], &vertices[1]),
        builder::line(&vertices[1], &vertices[2]),
        builder::line(&vertices[2], &vertices[3]),
    ]
    .into();
    let shell: Shell = builder::revolve_wire(
        &profile,
        Point3::origin(),
        Vector3::unit_z(),
        builder::SweepAngle::Closed,
        2,
    );
    let host = Solid::try_new(vec![shell])?;
    anyhow::ensure!(
        host.is_geometric_consistent(),
        "axis-aware revolved host is not geometrically consistent"
    );

    let cutter: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.0, -0.25, 1.0),
        Point3::new(2.0001, 0.25, 2.0001),
    ]));

    let result = monstertruck_solid::difference(&host, &cutter, TOL)?;
    anyhow::ensure!(
        !result.boundaries().is_empty(),
        "radial-slot difference returned an empty solid"
    );
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "radial-slot difference returned inconsistent topology"
    );
    Ok(())
}
