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
fn exact_line_arc_line_face_sweep_is_one_closed_solid() -> anyhow::Result<()> {
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
    assert_eq!(solid.boundaries()[0].len(), 14);
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
    assert_eq!(solid.boundaries()[0].len(), 16);
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
