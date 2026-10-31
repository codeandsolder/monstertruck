use super::common::*;
use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::*;

#[test]
fn step_redox_122_two_core_union_regression() -> Result<()> {
    fn frame_matrix(origin: Point3, x: Vector3, y: Vector3, z: Vector3) -> Matrix4 {
        Matrix4::from_cols(
            Vector4::new(x.x, x.y, x.z, 0.0),
            Vector4::new(y.x, y.y, y.z, 0.0),
            Vector4::new(z.x, z.y, z.z, 0.0),
            Vector4::new(origin.x, origin.y, origin.z, 1.0),
        )
    }

    let small_uv = [
        [0.0, 0.0],
        [0.6000000000000001, 0.0],
        [0.6000000000000001, 0.3966727112955528],
        [0.399999999999999, 0.5966727112955539],
        [0.20000000000000018, 0.5966727112955539],
        [0.0, 0.3966727112955546],
    ];
    let small_vertices = small_uv.map(|[x, y]| builder::vertex(Point3::new(x, y, 0.0)));
    let small_wire: Wire = (0..small_vertices.len())
        .map(|i| {
            builder::line(
                &small_vertices[i],
                &small_vertices[(i + 1) % small_vertices.len()],
            )
        })
        .collect();
    let small = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![small_wire],
        Vector3::unit_z() * 0.4,
    )?;
    let small = builder::transformed(
        &small,
        frame_matrix(
            Point3::new(3.5099999999999993, 3.9872458642577024, -0.12608142444673545),
            Vector3::new(1.0, -3.7854716527663227e-31, 1.0842021724855044e-16),
            Vector3::new(0.0, 1.0, 3.4914813388431334e-15),
            Vector3::new(1.0842021724855044e-16, -3.4914813388431334e-15, 1.0),
        ),
    );

    let main_uv = [
        [0.0, -0.40000000000000013],
        [1.3911086970663307, -0.3999999999999525],
        [2.0849019179718615, -0.6798257387093453],
        [2.8118807420625336, -1.3801742612906376],
        [3.5056739629680655, -1.659999999999988],
        [3.90234667426362, -1.6599999999999866],
        [3.9023466742636193, -1.2599999999999865],
        [3.5056739629680638, -1.259999999999988],
        [3.089398030424745, -1.0921045567743777],
        [2.3624192063340725, -0.39175603419308525],
        [1.3911086970663293, 0.0],
        [0.0, 0.0],
    ];
    let main_vertices = main_uv.map(|[x, y]| builder::vertex(Point3::new(x, y, 0.0)));
    let arc_mid = |center: [f64; 2], radius: f64, start: f64, end: f64| {
        let angle = (start + end) * 0.5;
        Point3::new(
            center[0] + radius * angle.cos(),
            center[1] + radius * angle.sin(),
            0.0,
        )
    };
    let main_wire: Wire = vec![
        builder::line(&main_vertices[0], &main_vertices[1]),
        builder::circle_arc::<Curve>(
            &main_vertices[1],
            &main_vertices[2],
            arc_mid(
                [1.391108697066334, -1.39999999999995],
                0.9999999999999998,
                1.5707963267948999,
                0.8040534579990185,
            ),
        )?,
        builder::line(&main_vertices[2], &main_vertices[3]),
        builder::circle_arc::<Curve>(
            &main_vertices[3],
            &main_vertices[4],
            arc_mid(
                [3.505673962968062, -0.659999999999879],
                0.9999999999999998,
                -2.3375391955907756,
                -1.570796326794893,
            ),
        )?,
        builder::line(&main_vertices[4], &main_vertices[5]),
        builder::line(&main_vertices[5], &main_vertices[6]),
        builder::line(&main_vertices[6], &main_vertices[7]),
        builder::circle_arc::<Curve>(
            &main_vertices[7],
            &main_vertices[8],
            arc_mid(
                [3.505673962968062, -0.65999999999879],
                0.599999999999999,
                -1.5707963267948937,
                -2.337539195590775,
            ),
        )?,
        builder::line(&main_vertices[8], &main_vertices[9]),
        builder::circle_arc::<Curve>(
            &main_vertices[9],
            &main_vertices[10],
            arc_mid(
                [1.391108697066334, -1.39999999999995],
                1.399999999999997,
                0.8040534579990186,
                1.5707963267948999,
            ),
        )?,
        builder::line(&main_vertices[10], &main_vertices[11]),
        builder::line(&main_vertices[11], &main_vertices[0]),
    ]
    .into();
    let main = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![main_wire],
        Vector3::unit_z() * 0.6000000000000001,
    )?;
    let main = builder::transformed(
        &main,
        frame_matrix(
            Point3::new(3.5099999999999993, 0.48157190128963684, 1.5339185755532525),
            Vector3::new(3.7854716527663227e-31, 1.0, 4.10421658981154e-47),
            Vector3::new(-1.0842021724855044e-16, 0.0, 1.0),
            Vector3::new(1.0, -3.7854716527663227e-31, 1.0842021724855044e-16),
        ),
    );

    let result = monstertruck_solid::or(
        &main,
        &small,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(result.is_geometric_consistent(), "union must be consistent");
    let stats = assert_union_matches_operand_occupancy(
        "step-redox case 122 core union",
        &main,
        &small,
        &result,
    )?;
    anyhow::ensure!(
        stats.lhs_only > 0 && stats.rhs_only > 0,
        "case 122 oracle must sample material contributed uniquely by both operands: {stats:?}"
    );
    Ok(())
}

#[test]
fn step_redox_38_three_core_union_regression() -> Result<()> {
    fn frame_matrix(origin: Point3, x: Vector3, y: Vector3, z: Vector3) -> Matrix4 {
        Matrix4::from_cols(
            Vector4::new(x.x, x.y, x.z, 0.0),
            Vector4::new(y.x, y.y, y.z, 0.0),
            Vector4::new(z.x, z.y, z.z, 0.0),
            Vector4::new(origin.x, origin.y, origin.z, 1.0),
        )
    }
    let point_on_arc = |center: [f64; 2], radius: f64, angle: f64| {
        Point3::new(
            center[0] + radius * angle.cos(),
            center[1] + radius * angle.sin(),
            0.0,
        )
    };
    let arc_mid = |center: [f64; 2], radius: f64, start: f64, end: f64| {
        point_on_arc(center, radius, (start + end) * 0.5)
    };

    let plate_starts = [
        Point3::new(0.0, 0.0, 0.0),
        point_on_arc(
            [0.5, -0.19999999999999996],
            0.2000000000000001,
            std::f64::consts::FRAC_PI_2,
        ),
        Point3::new(0.7, -0.19999999999999996, 0.0),
        Point3::new(0.7, -1.7, 0.0),
        Point3::new(2.7, -1.7000000000000002, 0.0),
        point_on_arc(
            [2.9, -0.19999999999999996],
            0.2000000000000001,
            std::f64::consts::PI,
        ),
        Point3::new(2.8999999999999995, 0.0, 0.0),
        Point3::new(3.4, 0.0, 0.0),
        Point3::new(3.4, 2.5, 0.0),
        Point3::new(0.0, 2.5, 0.0),
    ];
    let plate_vertices = plate_starts.map(builder::vertex);
    let plate_wire: Wire = vec![
        builder::line(&plate_vertices[0], &plate_vertices[1]),
        builder::circle_arc::<Curve>(
            &plate_vertices[1],
            &plate_vertices[2],
            arc_mid(
                [0.5, -0.19999999999999996],
                0.2000000000000001,
                std::f64::consts::FRAC_PI_2,
                0.0,
            ),
        )?,
        builder::line(&plate_vertices[2], &plate_vertices[3]),
        builder::line(&plate_vertices[3], &plate_vertices[4]),
        builder::line(&plate_vertices[4], &plate_vertices[5]),
        builder::circle_arc::<Curve>(
            &plate_vertices[5],
            &plate_vertices[6],
            arc_mid(
                [2.9, -0.19999999999999996],
                0.2000000000000001,
                std::f64::consts::PI,
                f64::from_bits(0x3ff9_21fb_5444_2d22),
            ),
        )?,
        builder::line(&plate_vertices[6], &plate_vertices[7]),
        builder::line(&plate_vertices[7], &plate_vertices[8]),
        builder::line(&plate_vertices[8], &plate_vertices[9]),
        builder::line(&plate_vertices[9], &plate_vertices[0]),
    ]
    .into();
    let plate = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![plate_wire],
        Vector3::unit_z() * 0.19999999999999973,
    )?;
    let plate = builder::transformed(
        &plate,
        frame_matrix(
            Point3::new(3.075000000000002, -1.3499999999999992, 2.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ),
    );

    let hook_starts = [
        point_on_arc(
            [0.30000000000000004, 0.0],
            0.3000000000000329,
            -std::f64::consts::PI,
        ),
        Point3::new(0.30000000000000004, -0.30000000000000027, 0.0),
        Point3::new(1.8, -0.30000000000000027, 0.0),
        Point3::new(1.8, -0.10000000000000053, 0.0),
        point_on_arc(
            [0.30000000000000004, 0.0],
            0.10000000000003278,
            4.712388980384706,
        ),
        Point3::new(0.2000000000000001, 0.0, 0.0),
        Point3::new(0.2000000000000001, 0.39999999999999813, 0.0),
        Point3::new(0.0, 0.3999999999999977, 0.0),
    ];
    let hook_vertices = hook_starts.map(builder::vertex);
    let hook_wire: Wire = vec![
        builder::circle_arc::<Curve>(
            &hook_vertices[0],
            &hook_vertices[1],
            arc_mid(
                [0.30000000000000004, 0.0],
                0.3000000000000329,
                -std::f64::consts::PI,
                -std::f64::consts::FRAC_PI_2,
            ),
        )?,
        builder::line(&hook_vertices[1], &hook_vertices[2]),
        builder::line(&hook_vertices[2], &hook_vertices[3]),
        builder::line(&hook_vertices[3], &hook_vertices[4]),
        builder::circle_arc::<Curve>(
            &hook_vertices[4],
            &hook_vertices[5],
            arc_mid(
                [0.30000000000000004, 0.0],
                0.10000000000003278,
                4.712388980384706,
                std::f64::consts::PI,
            ),
        )?,
        builder::line(&hook_vertices[5], &hook_vertices[6]),
        builder::line(&hook_vertices[6], &hook_vertices[7]),
        builder::line(&hook_vertices[7], &hook_vertices[0]),
    ]
    .into();
    let hook = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![hook_wire],
        Vector3::unit_z() * 2.0,
    )?;
    let hook = builder::transformed(
        &hook,
        frame_matrix(
            Point3::new(3.375000000000002, -0.6499999999999994, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, -1.2046690805394516e-16, 0.0),
            Vector3::new(1.2046690805394516e-16, 1.0, 0.0),
        ),
    );

    let two_core = monstertruck_solid::or(
        &plate,
        &hook,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        two_core.is_geometric_consistent(),
        "two-core union must be consistent"
    );
    let two_core_stats = assert_union_matches_operand_occupancy(
        "step-redox case 38 two-core union",
        &plate,
        &hook,
        &two_core,
    )?;
    anyhow::ensure!(
        two_core_stats.lhs_only > 0 && two_core_stats.rhs_only > 0,
        "case 38 two-core oracle must see unique material from both operands: {two_core_stats:?}"
    );

    // Missing third strict core from the actual Step Redox admission. The historical
    // two-core regression above exercises B+C, but the recovered solid is A+B+C.
    // A and B are disjoint; `or(A, B)` therefore has two boundary shells, and
    // `or(or(A, B), C)` must join both through C without misclassifying either
    // coplanar interface.
    let rounded_starts = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.3999999999999977, 0.0, 0.0),
        Point3::new(0.5999999999999979, 0.2, 0.0),
        Point3::new(0.5999999999999979, 1.7999999999999998, 0.0),
        Point3::new(0.3999999999999977, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    let rounded_vertices = rounded_starts.map(builder::vertex);
    let rounded_wire: Wire = vec![
        builder::line(&rounded_vertices[0], &rounded_vertices[1]),
        builder::circle_arc::<Curve>(
            &rounded_vertices[1],
            &rounded_vertices[2],
            arc_mid(
                [0.39999999999999813, 0.2],
                0.2000000000000001,
                f64::from_bits(0xbff9_21fb_5444_2d22),
                0.0,
            ),
        )?,
        builder::line(&rounded_vertices[2], &rounded_vertices[3]),
        builder::circle_arc::<Curve>(
            &rounded_vertices[3],
            &rounded_vertices[4],
            arc_mid(
                [0.39999999999999813, 1.7999999999999998],
                0.2000000000000001,
                0.0,
                f64::from_bits(0x3ff9_21fb_5444_2d22),
            ),
        )?,
        builder::line(&rounded_vertices[4], &rounded_vertices[5]),
        builder::line(&rounded_vertices[5], &rounded_vertices[0]),
    ]
    .into();
    let rounded = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![rounded_wire],
        Vector3::unit_z() * 0.2000000000000001,
    )?;
    let rounded = builder::transformed(
        &rounded,
        frame_matrix(
            Point3::new(3.375000000000002, -0.6499999999999994, 0.0),
            Vector3::unit_x(),
            Vector3::unit_y(),
            Vector3::unit_z(),
        ),
    );

    let disconnected = monstertruck_solid::or(
        &rounded,
        &plate,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        disconnected.is_geometric_consistent(),
        "A+B disconnected union must be consistent"
    );
    assert_union_matches_operand_occupancy(
        "step-redox case 38 disconnected A+B union",
        &rounded,
        &plate,
        &disconnected,
    )?;
    anyhow::ensure!(
        disconnected.boundaries().len() == 2,
        "A+B must remain two disconnected boundary shells"
    );

    let result = monstertruck_solid::or(
        &disconnected,
        &hook,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "three-core union must be consistent"
    );
    let three_core_stats = assert_union_matches_operand_occupancy(
        "step-redox case 38 three-core union",
        &disconnected,
        &hook,
        &result,
    )?;
    anyhow::ensure!(
        three_core_stats.lhs_only > 0 && three_core_stats.rhs_only > 0,
        "case 38 three-core oracle must see retained A+B material and hook-only material: {three_core_stats:?}"
    );
    Ok(())
}

#[test]
fn step_redox_173_edge_touching_cores_remain_distinct() -> Result<()> {
    fn frame_matrix(origin: Point3, x: Vector3, y: Vector3, z: Vector3) -> Matrix4 {
        Matrix4::from_cols(
            Vector4::new(x.x, x.y, x.z, 0.0),
            Vector4::new(y.x, y.y, y.z, 0.0),
            Vector4::new(z.x, z.y, z.z, 0.0),
            Vector4::new(origin.x, origin.y, origin.z, 1.0),
        )
    }

    let plate_uv = [
        [0.0, -2.8000000000000003],
        [0.5499999999999994, -2.8000000000000003],
        [0.5499999999999994, -2.999999999999996],
        [2.0499999999999994, -2.999999999999996],
        [2.0499999999999994, -2.8000000000000003],
        [2.6000000000000005, -2.8000000000000007],
        [2.6, 0.0],
        [0.0, 0.0],
    ];
    let plate_vertices = plate_uv.map(|[x, y]| builder::vertex(Point3::new(x, y, 0.0)));
    let plate_wire: Wire = (0..plate_vertices.len())
        .map(|index| {
            builder::line(
                &plate_vertices[index],
                &plate_vertices[(index + 1) % plate_vertices.len()],
            )
        })
        .collect();
    let plate = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![plate_wire],
        Vector3::unit_z() * 0.20000000000000018,
    )?;
    let plate = builder::transformed(
        &plate,
        frame_matrix(
            Point3::new(5.949999999999996, -1.9849999999999999, 3.099999999999997),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(5.782411586589357e-16, 0.0, 1.0),
            Vector3::new(1.0, 0.0, -5.782411586589357e-16),
        ),
    );

    let rounded_points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.10000000000000009, -0.10000000000000053, 0.0),
        Point3::new(1.4, -0.10000000000000053, 0.0),
        Point3::new(1.5, 0.0, 0.0),
        Point3::new(1.5, 0.510000000000006, 0.0),
        Point3::new(0.0, 0.510000000000006, 0.0),
    ];
    let rounded_vertices = rounded_points.map(builder::vertex);
    let rounded_wire: Wire = vec![
        builder::circle_arc::<Curve>(
            &rounded_vertices[0],
            &rounded_vertices[1],
            Point3::new(
                0.10000000000000009 - 0.10000000000000026 / f64::sqrt(2.0),
                -0.10000000000000026 / f64::sqrt(2.0),
                0.0,
            ),
        )?,
        builder::line(&rounded_vertices[1], &rounded_vertices[2]),
        builder::circle_arc::<Curve>(
            &rounded_vertices[2],
            &rounded_vertices[3],
            Point3::new(
                1.4 + 0.10000000000000026 / f64::sqrt(2.0),
                -0.10000000000000026 / f64::sqrt(2.0),
                0.0,
            ),
        )?,
        builder::line(&rounded_vertices[3], &rounded_vertices[4]),
        builder::line(&rounded_vertices[4], &rounded_vertices[5]),
        builder::line(&rounded_vertices[5], &rounded_vertices[0]),
    ]
    .into();
    let rounded = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![rounded_wire],
        Vector3::unit_z() * 0.2,
    )?;
    let rounded = builder::transformed(
        &rounded,
        frame_matrix(
            Point3::new(6.660000000000002, -1.4350000000000005, -0.09999999999999899),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 2.2217257632899447e-17),
            Vector3::new(2.2217257632899447e-17, 0.0, 1.0),
        ),
    );

    let result = monstertruck_solid::or(
        &plate,
        &rounded,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        result.boundaries().len() == 2,
        "edge-only contact must not weld two solids into one non-manifold shell"
    );
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "edge-touching union must remain geometrically consistent"
    );
    let stats = assert_union_matches_operand_occupancy(
        "step-redox case 173 edge-touching union",
        &plate,
        &rounded,
        &result,
    )?;
    anyhow::ensure!(
        stats.lhs_only > 0 && stats.rhs_only > 0,
        "case 173 oracle must see unique material from both edge-touching operands: {stats:?}"
    );
    Ok(())
}

#[test]
fn step_redox_97_nine_core_union_regression() -> Result<()> {
    fn frame_matrix(origin: Point3, x: Vector3, y: Vector3, z: Vector3) -> Matrix4 {
        Matrix4::from_cols(
            Vector4::new(x.x, x.y, x.z, 0.0),
            Vector4::new(y.x, y.y, y.z, 0.0),
            Vector4::new(z.x, z.y, z.z, 0.0),
            Vector4::new(origin.x, origin.y, origin.z, 1.0),
        )
    }
    let c1_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 5.590833366263048, 0.0)),
        builder::vertex(Point3::new(0.0, 5.590833366263048, 0.0)),
    ];
    let c1_wire: Wire = vec![
        builder::line(&c1_v[0], &c1_v[1]),
        builder::line(&c1_v[1], &c1_v[2]),
        builder::line(&c1_v[2], &c1_v[3]),
        builder::line(&c1_v[3], &c1_v[0]),
    ]
    .into();
    let c1 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c1_wire],
        Vector3::new(0.0, 0.0, 0.5000000000000009),
    )?;
    let c1 = builder::transformed(
        &c1,
        frame_matrix(
            Point3::new(-0.2500000000000011, -6.363043686338096, 3.5437360471411976),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, -1.3571966100182978e-16),
            Vector3::new(0.0, 1.3571966100182978e-16, 1.0),
        ),
    );
    let c2_v = [
        builder::vertex(Point3::new(
            -1.1102230246251565e-16,
            -1.1102230246251565e-16,
            0.0,
        )),
        builder::vertex(Point3::new(0.4979138337708573, -0.17383563328818652, 0.0)),
        builder::vertex(Point3::new(0.4979138337708573, 0.3261643667118138, 0.0)),
        builder::vertex(Point3::new(0.31119614610678603, 0.3913527291948844, 0.0)),
    ];
    let c2_wire: Wire = vec![
        builder::circle_arc::<Curve>(
            &c2_v[0],
            &c2_v[1],
            Point3::new(0.23422035012881803, -0.1291274621098183, 0.0),
        )?,
        builder::line(&c2_v[1], &c2_v[2]),
        builder::circle_arc::<Curve>(
            &c2_v[2],
            &c2_v[3],
            Point3::new(0.39902877740509246, 0.34292993090370183, 0.0),
        )?,
        builder::line(&c2_v[3], &c2_v[0]),
    ]
    .into();
    let c2 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c2_wire],
        Vector3::new(0.0, 0.0, 3.000000000000001),
    )?;
    let c2 = builder::transformed(
        &c2,
        frame_matrix(
            Point3::new(-0.2500000000000011, -5.631030843112577, 1.4756299721317527),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ),
    );
    let c3_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 4.086712518436875, 0.0)),
        builder::vertex(Point3::new(0.0, 4.086712518436875, 0.0)),
    ];
    let c3_wire: Wire = vec![
        builder::line(&c3_v[0], &c3_v[1]),
        builder::line(&c3_v[1], &c3_v[2]),
        builder::line(&c3_v[2], &c3_v[3]),
        builder::line(&c3_v[3], &c3_v[0]),
    ]
    .into();
    let c3 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c3_wire],
        Vector3::new(0.0, 0.0, 0.5000000000000004),
    )?;
    let c3 = builder::transformed(
        &c3,
        frame_matrix(
            Point3::new(-0.2500000000000011, -5.13311700934172, 1.3017943388435662),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, -3.9732422496187844e-16),
            Vector3::new(0.0, 3.9732422496187844e-16, 1.0),
        ),
    );
    let c4_v = [
        builder::vertex(Point3::new(
            -1.1102230246251565e-16,
            -9.797174393178825e-17,
            0.0,
        )),
        builder::vertex(Point3::new(0.3020861662291425, -0.6261643667118149, 0.0)),
        builder::vertex(Point3::new(0.613282312335928, -0.23481163751692993, 0.0)),
        builder::vertex(Point3::new(0.5, 0.0, 0.0)),
    ];
    let c4_wire: Wire = vec![
        builder::circle_arc::<Curve>(
            &c4_v[0],
            &c4_v[1],
            Point3::new(0.07946857562744503, -0.34761252349657545, 0.0),
        )?,
        builder::line(&c4_v[1], &c4_v[2]),
        builder::circle_arc::<Curve>(
            &c4_v[2],
            &c4_v[3],
            Point3::new(0.5298007158602916, -0.1303546963112155, 0.0),
        )?,
        builder::line(&c4_v[3], &c4_v[0]),
    ]
    .into();
    let c4 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c4_wire],
        Vector3::new(0.0, 0.0, 3.000000000000001),
    )?;
    let c4 = builder::transformed(
        &c4,
        frame_matrix(
            Point3::new(-0.2500000000000011, -7.163043686338096, 3.0798083229611994),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ),
    );
    let c5_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(0.5, 0.0, 0.0)),
        builder::vertex(Point3::new(0.5, 0.16392772417999835, 0.0)),
        builder::vertex(Point3::new(0.0, 0.16392772417999835, 0.0)),
    ];
    let c5_wire: Wire = vec![
        builder::line(&c5_v[0], &c5_v[1]),
        builder::line(&c5_v[1], &c5_v[2]),
        builder::line(&c5_v[2], &c5_v[3]),
        builder::line(&c5_v[3], &c5_v[0]),
    ]
    .into();
    let c5 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c5_wire],
        Vector3::new(0.0, 0.0, 3.000000000000001),
    )?;
    let c5 = builder::transformed(
        &c5,
        frame_matrix(
            Point3::new(-0.2500000000000011, -7.163043686338096, 3.0798083229611994),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ),
    );
    let c6_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(0.4999999999999999, 3.673940397442059e-17, 0.0)),
        builder::vertex(Point3::new(0.7999999999999998, 0.2999999999999998, 0.0)),
        builder::vertex(Point3::new(0.7999999999999998, 0.8000000000000004, 0.0)),
    ];
    let c6_wire: Wire = vec![
        builder::line(&c6_v[0], &c6_v[1]),
        builder::circle_arc::<Curve>(
            &c6_v[1],
            &c6_v[2],
            Point3::new(0.5878679656440357, 0.21213203435596423, 0.0),
        )?,
        builder::line(&c6_v[2], &c6_v[3]),
        builder::circle_arc::<Curve>(
            &c6_v[3],
            &c6_v[0],
            Point3::new(0.2343145750507616, 0.5656854249492383, 0.0),
        )?,
    ]
    .into();
    let c6 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c6_wire],
        Vector3::new(0.0, 0.0, 3.000000000000001),
    )?;
    let c6 = builder::transformed(
        &c6,
        frame_matrix(
            Point3::new(-0.2500000000000011, -7.163043686338096, 3.243736047141198),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ),
    );
    let c7_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(1.2299266769963761, -0.9780139841176319, 0.0)),
        builder::vertex(Point3::new(1.5411228231031622, -0.5866612549227475, 0.0)),
        builder::vertex(Point3::new(0.31119614610678514, 0.3913527291948844, 0.0)),
    ];
    let c7_wire: Wire = vec![
        builder::line(&c7_v[0], &c7_v[1]),
        builder::line(&c7_v[1], &c7_v[2]),
        builder::line(&c7_v[2], &c7_v[3]),
        builder::line(&c7_v[3], &c7_v[0]),
    ]
    .into();
    let c7 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c7_wire],
        Vector3::new(0.0, 0.0, 3.000000000000001),
    )?;
    let c7 = builder::transformed(
        &c7,
        frame_matrix(
            Point3::new(-0.2500000000000011, -6.860957520108953, 2.4536439562493846),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        ),
    );
    let c8_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 1.5713786888961527, 0.0)),
        builder::vertex(Point3::new(0.0, 1.5713786888961527, 0.0)),
    ];
    let c8_wire: Wire = vec![
        builder::line(&c8_v[0], &c8_v[1]),
        builder::line(&c8_v[1], &c8_v[2]),
        builder::line(&c8_v[2], &c8_v[3]),
        builder::line(&c8_v[3], &c8_v[0]),
    ]
    .into();
    let c8 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c8_wire],
        Vector3::new(0.0, 0.0, 0.5),
    )?;
    let c8 = builder::transformed(
        &c8,
        frame_matrix(
            Point3::new(-0.2500000000000011, -6.860957520108953, 2.4536439562493846),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.7827054583897682, -0.6223922922135708),
            Vector3::new(0.0, 0.6223922922135708, 0.7827054583897682),
        ),
    );
    let c9_v = [
        builder::vertex(Point3::new(0.0, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 0.0, 0.0)),
        builder::vertex(Point3::new(3.000000000000001, 0.16392772417999835, 0.0)),
        builder::vertex(Point3::new(0.0, 0.16392772417999835, 0.0)),
    ];
    let c9_wire: Wire = vec![
        builder::line(&c9_v[0], &c9_v[1]),
        builder::line(&c9_v[1], &c9_v[2]),
        builder::line(&c9_v[2], &c9_v[3]),
        builder::line(&c9_v[3], &c9_v[0]),
    ]
    .into();
    let c9 = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![c9_wire],
        Vector3::new(0.0, 0.0, 0.5),
    )?;
    let c9 = builder::transformed(
        &c9,
        frame_matrix(
            Point3::new(-0.2500000000000011, -7.163043686338096, 3.243736047141198),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, -1.0),
            Vector3::new(0.0, 1.0, 0.0),
        ),
    );
    let mut result = c1;
    for (index, rhs) in [&c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9]
        .into_iter()
        .enumerate()
    {
        let next = monstertruck_solid::or(
            &result,
            rhs,
            monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
        )
        .map_err(|error| anyhow::anyhow!("case 97 union child {} failed: {error}", index + 2))?;
        assert_union_matches_operand_occupancy(
            &format!("step-redox case 97 union child {}", index + 2),
            &result,
            rhs,
            &next,
        )?;
        result = next;
    }
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "case 97 nine-core union must be consistent"
    );
    anyhow::ensure!(
        result.boundaries().len() == 1,
        "case 97 union must produce one boundary shell"
    );
    Ok(())
}
