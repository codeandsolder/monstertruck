use super::common::*;
use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::*;

#[test]
fn union_of_overlapping_cubes() -> Result<()> {
    // Two unit cubes overlapping in a `0.5` cube: `2 - 0.5^3 = 1.875`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::or(&a, &b, monstertruck_solid::GeometryPolicy::standard(TOL))?;
    assert_solid("union of overlapping cubes", &result, 1.875)
}

#[test]
fn intersection_of_overlapping_cubes() -> Result<()> {
    // The shared `0.5` cube: `0.5^3 = 0.125`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result =
        monstertruck_solid::and(&a, &b, monstertruck_solid::GeometryPolicy::standard(TOL))?;
    assert_solid("intersection of overlapping cubes", &result, 0.125)
}

#[test]
fn difference_of_overlapping_cubes() -> Result<()> {
    // First cube minus the shared region: `1 - 0.125 = 0.875`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result =
        monstertruck_solid::difference(&a, &b, monstertruck_solid::GeometryPolicy::standard(TOL))?;
    assert_solid("difference of overlapping cubes", &result, 0.875)
}

#[test]
fn cube_minus_column() -> Result<()> {
    // Unit cube with a `0.4 x 0.4` column punched through: `1 - 0.4^2 = 0.84`.
    let cube = unit_cube(Point3::origin());
    let column = square_column(0.5, 0.5, 0.2);
    let result = monstertruck_solid::difference(
        &cube,
        &column,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        result.boundaries().len() == 1,
        "through-column subtraction must produce one connected boundary shell, got {}",
        result.boundaries().len(),
    );
    anyhow::ensure!(
        result.face_iter().count() == 10,
        "through-column subtraction must produce 10 faces, got {}",
        result.face_iter().count(),
    );
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
    let profile_host = profile_box(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 6.0, 2.0))?;
    let profile_cutter = profile_box(Point3::new(3.0, 2.0, -1.0), Point3::new(7.0, 4.0, 3.0))?;
    let profile_cutter_at_origin =
        profile_box(Point3::new(3.0, 2.0, 0.0), Point3::new(7.0, 4.0, 4.0))?;
    let transformed_profile_cutter = builder::transformed(
        &profile_cutter_at_origin,
        Matrix4::from_translation(Vector3::new(0.0, 0.0, -1.0)),
    );

    for (label, solid) in [
        ("primitive_host", &primitive_host),
        ("primitive_cutter", &primitive_cutter),
        ("profile_host", &profile_host),
        ("profile_cutter", &profile_cutter),
        ("transformed_profile_cutter", &transformed_profile_cutter),
    ] {
        eprintln!(
            "input_shell {label} condition={:?} consistent={}",
            solid.boundaries()[0].shell_condition(),
            solid.is_geometric_consistent(),
        );
    }

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
        match monstertruck_solid::difference(
            host,
            cutter,
            monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
        ) {
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
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
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
    let result = monstertruck_solid::difference(
        &outer,
        &inner,
        monstertruck_solid::GeometryPolicy::standard(TOL),
    )?;
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
    )?;
    let host = Solid::new(vec![shell])?;
    anyhow::ensure!(
        host.is_geometric_consistent(),
        "axis-aware revolved host is not geometrically consistent"
    );

    let cutter: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.0, -0.25, 1.0),
        Point3::new(2.0001, 0.25, 2.0001),
    ]));

    let result = monstertruck_solid::difference(
        &host,
        &cutter,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
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

#[test]
fn union_of_end_to_end_cylinders_with_full_circular_interface() -> Result<()> {
    let a = unit_cylinder_y(0.0)?;
    let b = unit_cylinder_y(1.0)?;
    let result = monstertruck_solid::or(
        &a,
        &b,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(result.boundaries().len() == 1, "union must have one shell");
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "union must be geometrically consistent"
    );
    Ok(())
}

#[test]
fn union_of_coaxial_overlapping_boxes_with_coplanar_sides() -> Result<()> {
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.0, 0.0, 0.5));
    let result =
        monstertruck_solid::or(&a, &b, monstertruck_solid::GeometryPolicy::standard(1.0e-6))?;
    assert_solid("union of coaxial overlapping boxes", &result, 1.5)
}

#[test]
fn intersection_of_coaxial_overlapping_boxes_with_coplanar_sides() -> Result<()> {
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.0, 0.0, 0.5));
    let result =
        monstertruck_solid::and(&a, &b, monstertruck_solid::GeometryPolicy::standard(1.0e-6))?;
    assert_solid("intersection of coaxial overlapping boxes", &result, 0.5)
}

/// Two volumes touch only across a partially overlapping planar interface.
/// Neither contact-face polygon contains a vertex of the other polygon:
/// their footprints cross. The interface must be imprinted and removed.
#[test]
fn union_of_orthogonal_boxes_touching_across_partial_face() -> Result<()> {
    let base: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 2.0, 1.0),
    ]));
    let cap: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.5, -1.0, 1.0),
        Point3::new(1.5, 3.0, 2.0),
    ]));
    for (left, right) in [(&base, &cap), (&cap, &base)] {
        let joined = monstertruck_solid::or(
            left,
            right,
            monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
        )?;
        anyhow::ensure!(
            joined.boundaries().len() == 1,
            "partial-face union must be connected"
        );
        assert_solid("partial-face planar-contact union", &joined, 8.0)?;
    }
    Ok(())
}

/// A hooked prism touches the underside of a hollow rectangular extrusion
/// through a planar strip opening into its central aperture. A label-only
/// Boolean result can be a perfectly closed copy of *just the hook*, so
/// topology validity alone is not evidence that the entire union survived.
#[test]
fn union_of_hollow_tube_and_underslung_hook_preserves_both_volumes() -> Result<()> {
    fn polygon(points: &[Point3]) -> Wire {
        let vertices = builder::vertices(points.iter().copied());
        (0..vertices.len())
            .map(|i| builder::line(&vertices[i], &vertices[(i + 1) % vertices.len()]))
            .collect::<Vec<_>>()
            .into()
    }
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| {
        polygon(&[
            Point3::new(x0, y0, 0.0),
            Point3::new(x1, y0, 0.0),
            Point3::new(x1, y1, 0.0),
            Point3::new(x0, y1, 0.0),
        ])
    };
    let tube = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![rect(0.0, 0.0, 1.36, 1.36), rect(0.2, 0.2, 1.16, 1.16)],
        Vector3::new(0.0, 0.0, 2.5),
    )?;
    // Cross-section in YZ, with its upper contact strip on Z=0.
    let yz = [
        (0.0, 0.0),
        (0.0, -0.219_021_120_132_902),
        (0.58, -0.439_021_120_132_9),
        (0.58, -0.52),
        (0.78, -0.52),
        (0.78, -0.300_978_879_867_095),
        (0.2, -0.080_978_879_867_097),
        (0.2, 0.0),
    ];
    let hook = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![polygon(&yz.map(|(y, z)| Point3::new(0.43, y, z)))],
        Vector3::new(0.5, 0.0, 0.0),
    )?;
    let section_area = yz
        .iter()
        .enumerate()
        .map(|(i, &(y0, z0))| {
            let (y1, z1) = yz[(i + 1) % yz.len()];
            y0 * z1 - y1 * z0
        })
        .sum::<f64>()
        .abs()
        * 0.5;
    let expected = (1.36_f64.powi(2) - 0.96_f64.powi(2)) * 2.5 + section_area * 0.5;
    for (lhs, rhs) in [(&tube, &hook), (&hook, &tube)] {
        let result = monstertruck_solid::or(
            lhs,
            rhs,
            monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
        )?;
        anyhow::ensure!(result.boundaries().len() == 1);
        anyhow::ensure!(result.is_geometric_consistent());
        assert_solid("hollow-tube-and-hook union", &result, expected)?;
    }
    Ok(())
}
