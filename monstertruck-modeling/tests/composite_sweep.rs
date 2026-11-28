use monstertruck_modeling::*;

fn rectangular_face() -> anyhow::Result<Face> {
    let vertices = builder::vertices([
        (1.75, -0.25, 4.0),
        (2.25, -0.25, 4.0),
        (2.25, 0.25, 4.0),
        (1.75, 0.25, 4.0),
    ]);
    let wire: Wire = vec![
        builder::line(&vertices[0], &vertices[1]),
        builder::line(&vertices[1], &vertices[2]),
        builder::line(&vertices[2], &vertices[3]),
        builder::line(&vertices[3], &vertices[0]),
    ]
    .into();
    Ok(builder::try_attach_plane(vec![wire])?)
}

#[test]
fn exact_line_arc_line_face_sweep_is_closed_and_consistent() -> anyhow::Result<()> {
    let face = rectangular_face()?;
    let solid = builder::composite_sweep(
        &face,
        [
            builder::CompositeSweepSegment::Translation(Vector3::new(0.0, 0.0, -4.0)),
            builder::CompositeSweepSegment::Rotation {
                origin: Point3::origin(),
                axis: Vector3::unit_y(),
                angle: Rad(std::f64::consts::FRAC_PI_2),
                division: 1,
            },
            builder::CompositeSweepSegment::Translation(Vector3::new(-4.0, 0.0, 0.0)),
        ],
    )?;

    assert_eq!(solid.boundaries().len(), 1);
    assert_eq!(
        solid.boundaries()[0].shell_condition(),
        shell::ShellCondition::Closed
    );
    assert!(solid.is_geometric_consistent());
    Ok(())
}

#[test]
fn composite_sweep_rejects_empty_and_degenerate_paths() -> anyhow::Result<()> {
    let face = rectangular_face()?;
    let empty = builder::composite_sweep(&face, []);
    assert_eq!(empty.unwrap_err(), errors::Error::EmptyCompositeSweep);

    let zero = builder::composite_sweep(
        &face,
        [builder::CompositeSweepSegment::Translation(Vector3::new(
            0.0, 0.0, 0.0,
        ))],
    );
    assert_eq!(
        zero.unwrap_err(),
        errors::Error::InvalidCompositeSweepSegment
    );
    Ok(())
}

#[test]
fn exact_closed_composite_sweep_has_no_duplicate_caps() -> anyhow::Result<()> {
    let face = rectangular_face()?;
    let quarter = std::f64::consts::FRAC_PI_2;
    let solid = builder::composite_closed_sweep(
        &face,
        [
            builder::CompositeSweepSegment::Rotation {
                origin: Point3::origin(),
                axis: Vector3::unit_y(),
                angle: Rad(quarter),
                division: 1,
            },
            builder::CompositeSweepSegment::Rotation {
                origin: Point3::origin(),
                axis: Vector3::unit_y(),
                angle: Rad(quarter),
                division: 1,
            },
            builder::CompositeSweepSegment::Rotation {
                origin: Point3::origin(),
                axis: Vector3::unit_y(),
                angle: Rad(quarter),
                division: 1,
            },
            builder::CompositeSweepSegment::Rotation {
                origin: Point3::origin(),
                axis: Vector3::unit_y(),
                angle: Rad(quarter),
                division: 1,
            },
        ],
    )?;

    assert_eq!(solid.boundaries().len(), 1);
    // A closed sweep has no start/end cap. Check that geometric fact directly;
    // the number of side faces is free to change if adjacent pieces are merged.
    let start_plane_z = 4.0;
    assert!(solid.boundaries()[0].iter().all(|face| {
        !face
            .vertex_iter()
            .all(|vertex| (vertex.point().z - start_plane_z).abs() < KERNEL_EPSILON)
    }));
    assert_eq!(
        solid.boundaries()[0].shell_condition(),
        shell::ShellCondition::Closed
    );
    assert!(solid.is_geometric_consistent());
    Ok(())
}

#[test]
fn closed_composite_sweep_rejects_an_open_path() -> anyhow::Result<()> {
    let face = rectangular_face()?;
    let result = builder::composite_closed_sweep(
        &face,
        [builder::CompositeSweepSegment::Translation(
            Vector3::unit_z(),
        )],
    );
    assert_eq!(result.unwrap_err(), errors::Error::CompositeSweepNotClosed);
    Ok(())
}

#[test]
fn closed_composite_sweep_normalizes_final_subdivision_before_closure() -> anyhow::Result<()> {
    use builder::CompositeSweepSegment::Rotation;
    let face = rectangular_face()?;
    let quarter = std::f64::consts::FRAC_PI_2;
    let segment = Rotation {
        origin: Point3::origin(),
        axis: Vector3::unit_y(),
        angle: Rad(quarter),
        division: 1,
    };
    let build = |sign: f64| {
        builder::composite_closed_sweep(
            &face,
            [
                segment,
                segment,
                segment,
                Rotation {
                    origin: Point3::origin(),
                    axis: Vector3::unit_y() * sign,
                    angle: Rad(quarter * sign),
                    division: 3,
                },
            ],
        )
    };
    let solid = build(1.0)?;
    let reversed_axis = build(-1.0)?;
    for result in [&solid, &reversed_axis] {
        assert_eq!(
            result.boundaries()[0].shell_condition(),
            shell::ShellCondition::Closed
        );
        assert!(result.is_geometric_consistent());
    }
    let positive = &solid.boundaries()[0];
    let negative = &reversed_axis.boundaries()[0];
    assert_eq!(positive.len(), negative.len());
    for (a, b) in positive.iter().zip(negative.iter()) {
        let a = a.oriented_surface();
        let b = b.oriented_surface();
        assert_near!(a.evaluate(0.5, 0.1), b.evaluate(0.5, 0.1));
        assert_near!(a.normal(0.5, 0.1), b.normal(0.5, 0.1));
    }
    Ok(())
}

#[test]
fn composite_sweep_rejects_invalid_rotation_before_building_any_segment() -> anyhow::Result<()> {
    let face = rectangular_face()?;
    let result = builder::composite_sweep(
        &face,
        [
            builder::CompositeSweepSegment::Translation(Vector3::unit_z()),
            builder::CompositeSweepSegment::Rotation {
                origin: Point3::origin(),
                axis: Vector3::zero(),
                angle: Rad(0.5),
                division: 4,
            },
        ],
    );
    assert_eq!(
        result.unwrap_err(),
        errors::Error::InvalidCompositeSweepSegment
    );
    Ok(())
}
