use super::common::*;
use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::*;

#[test]
fn union_with_thin_overlapping_bezier_protrusion_keeps_host() -> Result<()> {
    let host = unit_cube(Point3::origin());
    let feature = bezier_disk_extrusion(0.5, 0.5, 0.999_99, 0.2, 0.100_01)?;
    let result = monstertruck_solid::or(
        &host,
        &feature,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "union must stay consistent"
    );
    let mesh = result.triangulation(0.005).to_polygon();
    let volume = mesh.volume();
    anyhow::ensure!(
        volume > 1.0 && volume < 1.02,
        "union volume {volume:.9} dropped the host or added unrelated material"
    );
    Ok(())
}

#[test]
fn union_with_thin_holed_profile_boss() -> Result<()> {
    fn rect_wire(x0: f64, y0: f64, x1: f64, y1: f64) -> Wire {
        let v = builder::vertices([
            Point3::new(x0, y0, 0.0),
            Point3::new(x1, y0, 0.0),
            Point3::new(x1, y1, 0.0),
            Point3::new(x0, y1, 0.0),
        ]);
        vec![
            builder::line(&v[0], &v[1]),
            builder::line(&v[1], &v[2]),
            builder::line(&v[2], &v[3]),
            builder::line(&v[3], &v[0]),
        ]
        .into()
    }

    let host = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![rect_wire(-6.0, -6.0, 6.0, 6.0)],
        Vector3::new(0.0, 0.0, 0.7),
    )?;
    let boss = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![rect_wire(0.0, 0.0, 1.0, 1.0), rect_wire(0.3, 0.3, 0.7, 0.7)],
        Vector3::new(0.0, 0.0, 0.01001),
    )?;
    let boss = builder::transformed(
        &boss,
        Matrix4::from_translation(Vector3::new(4.0, 4.0, 0.69999)),
    );
    let result = monstertruck_solid::or(
        &host,
        &boss,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    // Independent closed form: host = 12*12*0.7 = 100.8. The boss profile
    // area is 1 - 0.4^2 = 0.84; only 0.01000 of its 0.01001 height lies
    // above the host, so the union adds exactly 0.0084.
    assert_solid("thin holed profile boss", &result, 100.8 + 0.84 * 0.01)
}

#[test]
fn union_with_step_redox_116_thin_holed_bezier_boss() -> Result<()> {
    let host_vertices_0 = builder::vertices([
        Point3::new(-20.32, 5.154083333333313, 3.88408333333333),
        Point3::new(-20.32, 7.694083333333322, 3.88408333333334),
        Point3::new(-20.32, 7.694083333333313, 6.424083333333341),
        Point3::new(-20.32, 7.262283333333313, 6.42408333333333),
        Point3::new(-20.32, 7.262283333333315, 5.916083333333339),
        Point3::new(-20.32, 5.585883333333315, 5.916083333333332),
        Point3::new(-20.32, 5.585883333333313, 6.42408333333333),
        Point3::new(-20.32, 5.154083333333313, 6.42408333333333),
    ]);
    let host_wire_0: Wire = vec![
        builder::line(&host_vertices_0[0], &host_vertices_0[1]),
        builder::line(&host_vertices_0[1], &host_vertices_0[2]),
        builder::line(&host_vertices_0[2], &host_vertices_0[3]),
        builder::line(&host_vertices_0[3], &host_vertices_0[4]),
        builder::line(&host_vertices_0[4], &host_vertices_0[5]),
        builder::line(&host_vertices_0[5], &host_vertices_0[6]),
        builder::line(&host_vertices_0[6], &host_vertices_0[7]),
        builder::line(&host_vertices_0[7], &host_vertices_0[0]),
    ]
    .into();
    let host = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![host_wire_0],
        Vector3::new(40.64, 0.0, 0.0),
    )?;
    let feature_vertices_0 = builder::vertices([
        Point3::new(-17.603072985038086, 5.14408333333332, 4.411745403277998),
        Point3::new(-17.582277826711188, 5.14408333333332, 4.3602773864189),
        Point3::new(-17.52630419221638, 5.14408333333332, 4.341215157954013),
        Point3::new(-17.452134794185948, 5.14408333333332, 4.3843651114805136),
        Point3::new(-17.451008389776735, 5.14408333333332, 4.3843651114805136),
        Point3::new(-17.451008389776735, 5.14408333333332, 4.346933826494174),
        Point3::new(-17.408551608193655, 5.14408333333332, 4.346933826494174),
        Point3::new(-17.408551608193655, 5.14408333333332, 4.502724220953719),
        Point3::new(-17.49606456614995, 5.14408333333332, 4.594569503563332),
        Point3::new(-17.57941849244226, 5.14408333333332, 4.572041415376723),
        Point3::new(-17.57941849244226, 5.14408333333332, 4.530797684694626),
        Point3::new(-17.500050304829244, 5.14408333333332, 4.560084199338842),
        Point3::new(-17.451008389776735, 5.14408333333332, 4.500817998109315),
        Point3::new(-17.522751686002874, 5.14408333333332, 4.4907670049181485),
    ]);
    let feature_wire_0: Wire = vec![
        builder::bezier(
            &feature_vertices_0[0],
            &feature_vertices_0[1],
            vec![
                Point3::new(-17.603072985038086, 5.14408333333332, 4.390083780020675),
                Point3::new(-17.596141265595787, 5.14408333333332, 4.372927774400205),
            ],
        ),
        builder::bezier(
            &feature_vertices_0[1],
            &feature_vertices_0[2],
            vec![
                Point3::new(-17.568414387827744, 5.14408333333332, 4.347626998437596),
                Point3::new(-17.549785391826337, 5.14408333333332, 4.341215157954013),
            ],
        ),
        builder::bezier(
            &feature_vertices_0[2],
            &feature_vertices_0[3],
            vec![
                Point3::new(-17.493638464345675, 5.14408333333332, 4.341215157954013),
                Point3::new(-17.468857567339484, 5.14408333333332, 4.355598475796185),
            ],
        ),
        builder::line(&feature_vertices_0[3], &feature_vertices_0[4]),
        builder::line(&feature_vertices_0[4], &feature_vertices_0[5]),
        builder::line(&feature_vertices_0[5], &feature_vertices_0[6]),
        builder::line(&feature_vertices_0[6], &feature_vertices_0[7]),
        builder::bezier(
            &feature_vertices_0[7],
            &feature_vertices_0[8],
            vec![
                Point3::new(-17.408551608193655, 5.14408333333332, 4.564069938018136),
                Point3::new(-17.437751476343788, 5.14408333333332, 4.594569503563332),
            ],
        ),
        builder::bezier(
            &feature_vertices_0[8],
            &feature_vertices_0[9],
            vec![
                Point3::new(-17.527603889611438, 5.14408333333332, 4.594569503563332),
                Point3::new(-17.555417413873574, 5.14408333333332, 4.587117905162304),
            ],
        ),
        builder::line(&feature_vertices_0[9], &feature_vertices_0[10]),
        builder::bezier(
            &feature_vertices_0[10],
            &feature_vertices_0[11],
            vec![
                Point3::new(-17.5555040603665, 5.14408333333332, 4.550206499133537),
                Point3::new(-17.52907687999352, 5.14408333333332, 4.560084199338842),
            ],
        ),
        builder::bezier(
            &feature_vertices_0[11],
            &feature_vertices_0[12],
            vec![
                Point3::new(-17.467384576958548, 5.14408333333332, 4.560084199338842),
                Point3::new(-17.451008389776735, 5.14408333333332, 4.54032879892822),
            ],
        ),
        builder::line(&feature_vertices_0[12], &feature_vertices_0[13]),
        builder::bezier(
            &feature_vertices_0[13],
            &feature_vertices_0[0],
            vec![
                Point3::new(-17.57629921869224, 5.14408333333332, 4.483315406517134),
                Point3::new(-17.603072985038086, 5.14408333333332, 4.45697487263708),
            ],
        ),
    ]
    .into();
    let feature_vertices_1 = builder::vertices([
        Point3::new(-17.560009678003368, 5.14408333333332, 4.415211262997424),
        Point3::new(-17.54874563391005, 5.14408333333332, 4.445884121528482),
        Point3::new(-17.507675196216123, 5.14408333333332, 4.460787318330511),
        Point3::new(-17.451008389776735, 5.14408333333332, 4.468412209717388),
        Point3::new(-17.451008389776735, 5.14408333333332, 4.445190949585058),
        Point3::new(-17.469204153312337, 5.14408333333332, 4.395802448560835),
        Point3::new(-17.515300087601847, 5.14408333333332, 4.376047048150216),
        Point3::new(-17.54770587599377, 5.14408333333332, 4.386964506270652),
    ]);
    let feature_wire_1: Wire = vec![
        builder::bezier(
            &feature_vertices_1[0],
            &feature_vertices_1[1],
            vec![
                Point3::new(-17.560009678003368, 5.14408333333332, 4.428381529938598),
                Point3::new(-17.556197232311074, 5.14408333333332, 4.4386058161156265),
            ],
        ),
        builder::bezier(
            &feature_vertices_1[1],
            &feature_vertices_1[2],
            vec![
                Point3::new(-17.541207389017263, 5.14408333333332, 4.453162426943634),
                Point3::new(-17.527517243118513, 5.14408333333332, 4.458014630552212),
            ],
        ),
        builder::line(&feature_vertices_1[2], &feature_vertices_1[3]),
        builder::line(&feature_vertices_1[3], &feature_vertices_1[4]),
        builder::bezier(
            &feature_vertices_1[4],
            &feature_vertices_1[5],
            vec![
                Point3::new(-17.451008389776735, 5.14408333333332, 4.425262256188587),
                Point3::new(-17.45707364428861, 5.14408333333332, 4.408799422513839),
            ],
        ),
        builder::bezier(
            &feature_vertices_1[5],
            &feature_vertices_1[6],
            vec![
                Point3::new(-17.48124801584316, 5.14408333333332, 4.382632181619647),
                Point3::new(-17.496671091601605, 5.14408333333332, 4.376047048150216),
            ],
        ),
        builder::bezier(
            &feature_vertices_1[6],
            &feature_vertices_1[7],
            vec![
                Point3::new(-17.52873029402181, 5.14408333333332, 4.376047048150216),
                Point3::new(-17.539561105649323, 5.14408333333332, 4.379686200857799),
            ],
        ),
        builder::bezier(
            &feature_vertices_1[7],
            &feature_vertices_1[0],
            vec![
                Point3::new(-17.555937292831135, 5.14408333333332, 4.394242811685818),
                Point3::new(-17.560009678003368, 5.14408333333332, 4.403600632933563),
            ],
        ),
    ]
    .into();
    let feature = profile::solid_from_planar_profile::<Curve, Surface>(
        vec![feature_wire_0, feature_wire_1],
        Vector3::new(0.0, 0.010010000000000675, 0.0),
    )?;
    let result = monstertruck_solid::or(
        &host,
        &feature,
        monstertruck_solid::GeometryPolicy::standard(STRICT_TOL),
    )?;
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "step-redox 116 root1+microboss union must stay consistent"
    );
    let stats = assert_union_matches_operand_occupancy(
        "step-redox 116 root1+microboss",
        &host,
        &feature,
        &result,
    )?;
    anyhow::ensure!(
        stats.lhs_only > 0 && stats.rhs_only > 0,
        "step-redox 116 oracle must sample host-only and boss-only material: {stats:?}"
    );
    Ok(())
}
