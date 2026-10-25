use super::common::*;
use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::*;

#[test]
#[allow(
    clippy::excessive_precision,
    reason = "regression fixture preserves recovered STEP coordinates verbatim"
)]
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
        ),
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
        ),
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
        ),
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
        ),
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

    let result = monstertruck_solid::or(&main, &small, STRICT_TOL)?;
    anyhow::ensure!(result.is_geometric_consistent(), "union must be consistent");
    assert_solid(
        "step-redox case 122 core union",
        &result,
        1.0713717537932466,
    )
}

#[test]
#[allow(
    clippy::approx_constant,
    clippy::excessive_precision,
    reason = "regression fixture preserves recovered STEP coordinates verbatim"
)]
fn step_redox_38_two_core_union_regression() -> Result<()> {
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
            1.5707963267948966,
        ),
        Point3::new(0.7, -0.19999999999999996, 0.0),
        Point3::new(0.7, -1.7, 0.0),
        Point3::new(2.7, -1.7000000000000002, 0.0),
        point_on_arc(
            [2.9, -0.19999999999999996],
            0.2000000000000001,
            3.141592653589793,
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
                1.5707963267948966,
                0.0,
            ),
        ),
        builder::line(&plate_vertices[2], &plate_vertices[3]),
        builder::line(&plate_vertices[3], &plate_vertices[4]),
        builder::line(&plate_vertices[4], &plate_vertices[5]),
        builder::circle_arc::<Curve>(
            &plate_vertices[5],
            &plate_vertices[6],
            arc_mid(
                [2.9, -0.19999999999999996],
                0.2000000000000001,
                3.141592653589793,
                1.5707963267948988,
            ),
        ),
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
            -3.141592653589793,
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
                -3.141592653589793,
                -1.5707963267948966,
            ),
        ),
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
                3.141592653589793,
            ),
        ),
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

    let result = monstertruck_solid::or(&plate, &hook, STRICT_TOL)?;
    anyhow::ensure!(result.is_geometric_consistent(), "union must be consistent");
    assert_solid("step-redox case 38 core union", &result, 2.6690973355002896)
}
