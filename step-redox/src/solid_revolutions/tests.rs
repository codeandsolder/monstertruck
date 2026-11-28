#![expect(
    clippy::unreadable_literal,
    reason = "legacy numeric fixtures preserve source-derived expected coordinates"
)]
#![expect(
    clippy::suboptimal_flops,
    reason = "fixture arithmetic mirrors the geometric construction directly"
)]
#![expect(
    clippy::too_many_lines,
    reason = "large fixture constructors remain linear for reviewability"
)]
#![expect(
    clippy::assert_is_empty,
    reason = "test assertions intentionally state emptiness expectations explicitly"
)]
use super::*;

fn test_linear_profile_segment(
    face_index: usize,
    face: &FaceInfo,
    context: &TopologyContext<'_>,
    source_tolerance_mm: f64,
) -> Option<(Segment2, f64)> {
    linear_profile_segment(
        face_index,
        face,
        ProfileSegmentContext::without_angular_trims(
            context,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            source_tolerance_mm,
        ),
    )
}

fn profile_has_line(curves: &[RecoveredProfileCurve], first: [f64; 2], second: [f64; 2]) -> bool {
    curves.iter().any(|curve| {
        let RecoveredProfileCurve::Line {
            start_mm, end_mm, ..
        } = curve
        else {
            return false;
        };
        (distance2(*start_mm, first) <= 1.0e-12 && distance2(*end_mm, second) <= 1.0e-12)
            || (distance2(*start_mm, second) <= 1.0e-12 && distance2(*end_mm, first) <= 1.0e-12)
    })
}

fn profile_has_arc(
    curves: &[RecoveredProfileCurve],
    center: [f64; 2],
    radius: f64,
    first: [f64; 2],
    second: [f64; 2],
) -> bool {
    curves.iter().any(|curve| {
        let RecoveredProfileCurve::CircleArc {
            center_mm,
            radius_mm,
            ..
        } = curve
        else {
            return false;
        };
        let Some(start) = curve.start_point() else {
            return false;
        };
        let Some(end) = curve.end_point() else {
            return false;
        };
        distance2(*center_mm, center) <= 1.0e-12
            && (*radius_mm - radius).abs() <= 1.0e-12
            && ((distance2(start, first) <= 1.0e-12 && distance2(end, second) <= 1.0e-12)
                || (distance2(start, second) <= 1.0e-12 && distance2(end, first) <= 1.0e-12))
    })
}

#[test]
fn reads_representation_length_uncertainty_with_hard_cap() -> anyhow::Result<()> {
    fn tolerance_from(text: &str) -> anyhow::Result<Option<f64>> {
        let exchange = step_p21::parser::parse(text)?;
        let entities = &exchange
            .data
            .first()
            .ok_or_else(|| anyhow::anyhow!("parsed exchange has no DATA section"))?
            .entities;
        let index = build_index(entities);
        let solid_id = entities
            .iter()
            .find(|entity| {
                simple_record(entity).is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
            })
            .map(entity_id)
            .ok_or_else(|| anyhow::anyhow!("fixture has no MANIFOLD_SOLID_BREP"))?;
        Ok(source_tolerance_by_representation_item(entities, &index)
            .get(&solid_id)
            .copied())
    }

    let source = std::str::from_utf8(include_bytes!(
        "../../validation/fixtures/native_conical_frustum.step"
    ))?;
    assert_eq!(tolerance_from(source)?, Some(GEOM_TOL_MM));

    let widened = source.replacen("LENGTH_MEASURE(1.E-07)", "LENGTH_MEASURE(5.E-06)", 1);
    assert!(
        (tolerance_from(&widened)?
            .ok_or_else(|| anyhow::anyhow!("widened uncertainty missing"))?
            - 5.0e-6)
            .abs()
            <= 1.0e-15
    );

    let excessive = source.replacen("LENGTH_MEASURE(1.E-07)", "LENGTH_MEASURE(2.E-05)", 1);
    assert_eq!(tolerance_from(&excessive)?, None);

    let centimetres = source.replacen(
        "LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.)",
        "LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.CENTI.,.METRE.)",
        1,
    );
    assert!(
        (tolerance_from(&centimetres)?
            .ok_or_else(|| anyhow::anyhow!("centimetre uncertainty missing"))?
            - 1.0e-6)
            .abs()
            <= 1.0e-15
    );

    let ambiguous = source.replacen(
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117))",
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117,#117))",
        1,
    );
    assert_eq!(tolerance_from(&ambiguous)?, None);
    Ok(())
}

fn test_circle_edge(
    edge_id: u64,
    start_mm: [f64; 3],
    end_mm: [f64; 3],
    center_mm: [f64; 3],
    normal: [f64; 3],
    radius_mm: f64,
) -> brep::OrientedEdgeUse {
    brep::OrientedEdgeUse {
        oriented_edge_id: edge_id + 10_000,
        edge_id,
        curve_id: edge_id + 20_000,
        curve_same_sense: true,
        parameter_forward: true,
        start_vertex: edge_id + 30_000,
        end_vertex: edge_id + 40_000,
        start_mm,
        end_mm,
        support: CurveSupport::Circle(brep::CircleSupport {
            center_mm,
            normal,
            x_direction: [1.0, 0.0, 0.0],
            radius_mm,
        }),
    }
}

fn test_line_edge(edge_id: u64, start_mm: [f64; 3], end_mm: [f64; 3]) -> brep::OrientedEdgeUse {
    let delta = sub(end_mm, start_mm);
    let length = norm(delta);
    assert!(length > 1.0e-12, "test edge must have nonzero length");
    let direction = mul(delta, 1.0 / length);
    brep::OrientedEdgeUse {
        oriented_edge_id: edge_id + 10_000,
        edge_id,
        curve_id: edge_id + 20_000,
        curve_same_sense: true,
        parameter_forward: true,
        start_vertex: edge_id + 30_000,
        end_vertex: edge_id + 40_000,
        start_mm,
        end_mm,
        support: CurveSupport::Line(brep::LineSupport {
            origin_mm: start_mm,
            direction,
        }),
    }
}

fn test_face(surface: SurfaceSupport, edges: Vec<brep::OrientedEdgeUse>) -> FaceInfo {
    FaceInfo {
        surface,
        loops: vec![brep::FaceLoop {
            bound_id: 1,
            loop_id: 2,
            outer: true,
            orientation: true,
            edges,
        }],
    }
}

#[test]
fn plane_profile_accepts_bounded_noisy_circle_vertices() -> anyhow::Result<()> {
    let make_plane = |noise: f64| {
        let positive = [1.0, 0.0, noise];
        let negative = [-1.0, 0.0, 0.0];
        test_face(
            SurfaceSupport::Plane(brep::PlaneSupport {
                origin_mm: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                max_residual_mm: 0.0,
            }),
            vec![
                test_circle_edge(
                    800,
                    positive,
                    negative,
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    1.0,
                ),
                test_circle_edge(
                    801,
                    negative,
                    positive,
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    1.0,
                ),
            ],
        )
    };

    let noisy = make_plane(0.5 * REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    let faces = vec![noisy.clone()];
    let edge_faces = edge_face_map(&faces);
    let entities = Vec::new();
    let index = HashMap::new();
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    let (segment, residual) =
        test_linear_profile_segment(0, &noisy, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("expected bounded profile segment"))?;
    assert_eq!(segment.a, [0.0, 0.0]);
    assert_eq!(segment.b, [1.0, 0.0]);
    assert!((residual - 0.5 * REVOLUTION_SOURCE_SUPPORT_TOL_MM).abs() <= 1.0e-12);

    let too_noisy = make_plane(2.0 * REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    let faces = vec![too_noisy.clone()];
    let edge_faces = edge_face_map(&faces);
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    assert!(
        test_linear_profile_segment(0, &too_noisy, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .is_none()
    );
    Ok(())
}

#[test]
fn plane_profile_detects_axis_crossing_inside_radial_edge() -> anyhow::Result<()> {
    let make_face = |line_start: [f64; 3], line_end: [f64; 3]| {
        test_face(
            SurfaceSupport::Plane(brep::PlaneSupport {
                origin_mm: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                max_residual_mm: 0.0,
            }),
            vec![
                test_circle_edge(
                    820,
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    1.0,
                ),
                test_line_edge(821, line_start, line_end),
            ],
        )
    };

    let crossing = make_face([0.0, -0.25, 0.0], [0.0, 0.25, 0.0]);
    let faces = vec![crossing.clone()];
    let edge_faces = edge_face_map(&faces);
    let entities = Vec::new();
    let index = HashMap::new();
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    let (segment, _) =
        test_linear_profile_segment(0, &crossing, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("axis-crossing radial edge should prove a disk"))?;
    assert_eq!(segment.a, [0.0, 0.0]);
    assert_eq!(segment.b, [1.0, 0.0]);

    let one_sided = make_face([0.0, 0.25, 0.0], [0.0, 0.5, 0.0]);
    let faces = vec![one_sided.clone()];
    let edge_faces = edge_face_map(&faces);
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    assert!(
        test_linear_profile_segment(0, &one_sided, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .is_none()
    );
    Ok(())
}

#[test]
fn cylinder_profile_accepts_bounded_noisy_trim_vertices() -> anyhow::Result<()> {
    let make_cylinder = |noise: f64| {
        let bottom_pos = [1.0 + noise, 0.0, 0.0];
        let bottom_neg = [-1.0, 0.0, 0.0];
        let top_pos = [1.0 + noise, 0.0, 1.0];
        let top_neg = [-1.0, 0.0, 1.0];
        test_face(
            SurfaceSupport::Cylinder(brep::CylinderSupport {
                axis_origin_mm: [0.0, 0.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                x_direction: [1.0, 0.0, 0.0],
                radius_mm: 1.0,
            }),
            vec![
                test_circle_edge(
                    900,
                    bottom_pos,
                    bottom_neg,
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    1.0,
                ),
                test_line_edge(901, bottom_neg, top_neg),
                test_circle_edge(902, top_neg, top_pos, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], 1.0),
                test_line_edge(903, top_pos, bottom_pos),
            ],
        )
    };

    let empty_edge_faces = HashMap::new();
    let empty_entities = Vec::new();
    let empty_index = HashMap::new();

    let noisy = make_cylinder(0.5 * REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    let faces = vec![noisy.clone()];
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &empty_edge_faces,
        entities: &empty_entities,
        index: &empty_index,
    };
    let (segment, residual) =
        test_linear_profile_segment(0, &noisy, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("expected bounded profile segment"))?;
    assert_eq!(segment.a, [1.0, 0.0]);
    assert_eq!(segment.b, [1.0, 1.0]);
    assert!((residual - 0.5 * REVOLUTION_SOURCE_SUPPORT_TOL_MM).abs() <= 1.0e-12);

    let too_noisy = make_cylinder(2.0 * REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    let faces = vec![too_noisy.clone()];
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &empty_edge_faces,
        entities: &empty_entities,
        index: &empty_index,
    };
    assert!(
        test_linear_profile_segment(0, &too_noisy, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .is_none()
    );
    Ok(())
}

#[test]
fn cone_profile_accepts_proven_apex_and_bounded_source_noise() -> anyhow::Result<()> {
    let apex = [0.0, 0.0, 0.0];
    let base_pos = [1.0, 0.0, 1.0];
    let base_neg = [-1.0, 0.0, 1.0];
    let mut apex_line_a = test_line_edge(1000, apex, base_pos);
    let mut apex_line_b = test_line_edge(1002, base_neg, apex);
    apex_line_a.start_vertex = 42_424;
    apex_line_b.end_vertex = 42_424;
    let apex_face = test_face(
        SurfaceSupport::Cone(brep::ConeSupport {
            reference_origin_mm: [0.0, 0.0, 1.0],
            axis: [0.0, 0.0, 1.0],
            x_direction: [1.0, 0.0, 0.0],
            reference_radius_mm: 1.0,
            semi_angle_rad: std::f64::consts::FRAC_PI_4,
        }),
        vec![
            apex_line_a,
            test_circle_edge(
                1001,
                base_pos,
                base_neg,
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                1.0,
            ),
            apex_line_b,
        ],
    );
    let empty_edge_faces = HashMap::new();
    let empty_entities = Vec::new();
    let empty_index = HashMap::new();
    let faces = vec![apex_face.clone()];
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &empty_edge_faces,
        entities: &empty_entities,
        index: &empty_index,
    };
    let (segment, residual) =
        test_linear_profile_segment(0, &apex_face, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("expected bounded profile segment"))?;
    assert_eq!(segment.a, [0.0, 0.0]);
    assert_eq!(segment.b, [1.0, 1.0]);
    assert!(residual <= 1.0e-12);

    let make_frustum = |reference_radius_mm: f64| {
        let bottom_pos = [1.0, 0.0, 0.0];
        let bottom_neg = [-1.0, 0.0, 0.0];
        let top_pos = [2.0, 0.0, 1.0];
        let top_neg = [-2.0, 0.0, 1.0];
        test_face(
            SurfaceSupport::Cone(brep::ConeSupport {
                reference_origin_mm: [0.0, 0.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                x_direction: [1.0, 0.0, 0.0],
                reference_radius_mm,
                semi_angle_rad: std::f64::consts::FRAC_PI_4,
            }),
            vec![
                test_circle_edge(
                    1010,
                    bottom_pos,
                    bottom_neg,
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    1.0,
                ),
                test_line_edge(1011, bottom_neg, top_neg),
                test_circle_edge(
                    1012,
                    top_neg,
                    top_pos,
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, 1.0],
                    2.0,
                ),
                test_line_edge(1013, top_pos, bottom_pos),
            ],
        )
    };

    let noisy = make_frustum(1.0 + 0.5 * REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    let faces = vec![noisy.clone()];
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &empty_edge_faces,
        entities: &empty_entities,
        index: &empty_index,
    };
    let (segment, residual) =
        test_linear_profile_segment(0, &noisy, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("expected bounded profile segment"))?;
    assert_eq!(segment.a, [1.0, 0.0]);
    assert_eq!(segment.b, [2.0, 1.0]);
    assert!((residual - 0.5 * REVOLUTION_SOURCE_SUPPORT_TOL_MM).abs() <= 1.0e-12);

    let too_noisy = make_frustum(1.0 + 2.0 * REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    let faces = vec![too_noisy.clone()];
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &empty_edge_faces,
        entities: &empty_entities,
        index: &empty_index,
    };
    assert!(
        test_linear_profile_segment(0, &too_noisy, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .is_none()
    );

    let declared_tolerance_mm = 1.0e-6;
    let declared_noisy = make_frustum(1.0 + 0.5 * declared_tolerance_mm);
    let faces = vec![declared_noisy.clone()];
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &empty_edge_faces,
        entities: &empty_entities,
        index: &empty_index,
    };
    assert!(
        test_linear_profile_segment(
            0,
            &declared_noisy,
            &context,
            REVOLUTION_SOURCE_SUPPORT_TOL_MM
        )
        .is_none()
    );
    assert!(
        test_linear_profile_segment(0, &declared_noisy, &context, declared_tolerance_mm).is_some()
    );
    Ok(())
}

#[test]
fn cone_profile_uses_neighbor_constraint_plus_repeated_trim_circle() -> anyhow::Result<()> {
    let bottom_pos = [1.5, 0.0, 0.5000005];
    let bottom_neg = [-1.5, 0.0, 0.5000005];
    let top_pos = [2.0000004, 0.0, 1.0000004];
    let top_neg = [-2.0000004, 0.0, 1.0000004];

    let bottom = test_circle_edge(
        1100,
        bottom_pos,
        bottom_neg,
        [0.0, 0.0, 0.5000005],
        [0.0, 0.0, 1.0],
        1.5,
    );
    let top = test_circle_edge(
        1102,
        top_neg,
        top_pos,
        [0.0, 0.0, 1.0000004],
        [0.0, 0.0, 1.0],
        2.0000004,
    );
    let cone_face = test_face(
        SurfaceSupport::Cone(brep::ConeSupport {
            reference_origin_mm: [0.0, 0.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            x_direction: [1.0, 0.0, 0.0],
            reference_radius_mm: 1.0,
            semi_angle_rad: std::f64::consts::FRAC_PI_4,
        }),
        vec![
            bottom.clone(),
            test_line_edge(1101, bottom_neg, top_neg),
            top.clone(),
            test_line_edge(1103, top_pos, bottom_pos),
        ],
    );
    let cylinder_face = test_face(
        SurfaceSupport::Cylinder(brep::CylinderSupport {
            axis_origin_mm: [0.0, 0.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            x_direction: [1.0, 0.0, 0.0],
            radius_mm: 1.5,
        }),
        vec![bottom],
    );
    let plane_face = test_face(
        SurfaceSupport::Plane(brep::PlaneSupport {
            origin_mm: [0.0, 0.0, 1.0],
            normal: [0.0, 0.0, 1.0],
            max_residual_mm: 0.0,
        }),
        vec![top],
    );
    let faces = vec![cone_face.clone(), cylinder_face, plane_face];
    let edge_faces = edge_face_map(&faces);
    let entities = Vec::new();
    let index = HashMap::new();
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };

    assert!(
        test_linear_profile_segment(0, &cone_face, &context, REVOLUTION_SOURCE_SUPPORT_TOL_MM)
            .is_none()
    );

    let (segment, residual) = test_linear_profile_segment(0, &cone_face, &context, 1.0e-5)
        .ok_or_else(|| anyhow::anyhow!("expected recovered test geometry"))?;

    // The cylinder constrains radius while preserving the repeated trim-circle
    // axial coordinate; the plane constrains axial position while preserving
    // the repeated trim-circle radius. The cone support is evidence, not the
    // sole source of truth when its own parameters disagree within STEP uncertainty.
    assert_eq!(segment.a, [1.5, 0.5000005]);
    assert_eq!(segment.b, [2.0000004, 1.0]);
    assert!(residual > REVOLUTION_SOURCE_SUPPORT_TOL_MM);
    assert!(residual < 1.0e-5);
    Ok(())
}

fn split_hemisphere_faces(cylinder_radius_mm: f64) -> Vec<FaceInfo> {
    let sphere = brep::SphereSupport {
        center_mm: [0.0, 0.0, 0.0],
        axis: [0.0, 0.0, 1.0],
        x_direction: [1.0, 0.0, 0.0],
        radius_mm: 1.0,
    };
    let cylinder = brep::CylinderSupport {
        axis_origin_mm: [0.0, 0.0, 0.0],
        axis: [0.0, 0.0, 1.0],
        x_direction: [1.0, 0.0, 0.0],
        radius_mm: cylinder_radius_mm,
    };
    let plane = brep::PlaneSupport {
        origin_mm: [0.0, 0.0, 2.0],
        normal: [0.0, 0.0, 1.0],
        max_residual_mm: 0.0,
    };

    let equator_a = test_circle_edge(
        10,
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        1.0,
    );
    let equator_b = test_circle_edge(
        11,
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        1.0,
    );
    let sphere_seam_a = test_circle_edge(
        12,
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.0,
    );
    let sphere_seam_b = test_circle_edge(
        13,
        [0.0, 0.0, -1.0],
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.0,
    );
    let disk_a = test_circle_edge(
        14,
        [1.0, 0.0, 2.0],
        [-1.0, 0.0, 2.0],
        [0.0, 0.0, 2.0],
        [0.0, 0.0, 1.0],
        1.0,
    );
    let disk_b = test_circle_edge(
        15,
        [-1.0, 0.0, 2.0],
        [1.0, 0.0, 2.0],
        [0.0, 0.0, 2.0],
        [0.0, 0.0, 1.0],
        1.0,
    );
    let cylinder_seam_a = test_line_edge(16, [1.0, 0.0, 0.0], [1.0, 0.0, 2.0]);
    let cylinder_seam_b = test_line_edge(17, [-1.0, 0.0, 0.0], [-1.0, 0.0, 2.0]);

    vec![
        test_face(
            SurfaceSupport::Sphere(sphere),
            vec![
                equator_a.clone(),
                sphere_seam_a.clone(),
                sphere_seam_b.clone(),
            ],
        ),
        test_face(
            SurfaceSupport::Sphere(sphere),
            vec![equator_b.clone(), sphere_seam_a, sphere_seam_b],
        ),
        test_face(
            SurfaceSupport::Cylinder(cylinder),
            vec![
                equator_a,
                disk_a.clone(),
                cylinder_seam_a.clone(),
                cylinder_seam_b.clone(),
            ],
        ),
        test_face(
            SurfaceSupport::Cylinder(cylinder),
            vec![equator_b, disk_b.clone(), cylinder_seam_a, cylinder_seam_b],
        ),
        test_face(SurfaceSupport::Plane(plane), vec![disk_a, disk_b]),
    ]
}

fn edge_face_map(faces: &[FaceInfo]) -> HashMap<u64, Vec<usize>> {
    let mut edge_faces = HashMap::<u64, Vec<usize>>::new();
    for (face_index, face) in faces.iter().enumerate() {
        for edge in face.loops.iter().flat_map(|loop_| &loop_.edges) {
            edge_faces.entry(edge.edge_id).or_default().push(face_index);
        }
    }
    edge_faces
}

#[test]
fn recovers_split_hemispherical_end_topology() -> anyhow::Result<()> {
    let faces = split_hemisphere_faces(1.0);
    let edge_faces = edge_face_map(&faces);
    let entities = Vec::new();
    let index = HashMap::new();
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    let recovered = detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5], &faces, &context)
        .ok_or_else(|| anyhow::anyhow!("expected hemispherical recovery"))?;

    assert_eq!(recovered.profile_curves.len(), 4);
    assert!(recovered.max_residual_mm <= 1.0e-12);
    assert_eq!(recovered.axis_direction, [0.0, 0.0, 1.0]);

    assert!(profile_has_line(
        &recovered.profile_curves,
        [0.0, 2.0],
        [1.0, 2.0]
    ));
    assert!(profile_has_line(
        &recovered.profile_curves,
        [1.0, 2.0],
        [1.0, 0.0]
    ));
    assert!(profile_has_line(
        &recovered.profile_curves,
        [0.0, -1.0],
        [0.0, 2.0]
    ));
    assert!(profile_has_arc(
        &recovered.profile_curves,
        [0.0, 0.0],
        1.0,
        [1.0, 0.0],
        [0.0, -1.0],
    ));

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered)?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
        step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
    }

    let tampered = split_hemisphere_faces(1.01);
    let edge_faces = edge_face_map(&tampered);
    let context = TopologyContext {
        faces: &tampered,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    assert!(detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5], &tampered, &context).is_none());
    Ok(())
}

#[test]
fn mixed_curved_graph_recovers_split_sphere_arc() -> anyhow::Result<()> {
    let faces = split_hemisphere_faces(1.0);
    let edge_faces = edge_face_map(&faces);
    let entities = Vec::new();
    let index = HashMap::new();
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    let recovered = detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5], &faces, &context)
        .ok_or_else(|| anyhow::anyhow!("expected mixed-curved recovery"))?;
    assert_eq!(recovered.profile_curves.len(), 4);
    assert!(recovered.profile_curves.iter().any(|curve| matches!(
        curve,
        RecoveredProfileCurve::CircleArc {
            center_mm: [0.0, 0.0],
            radius_mm: 1.0,
            start_angle_rad,
            end_angle_rad,
            ..
        } if (start_angle_rad.abs() <= 1.0e-12
            && (*end_angle_rad + std::f64::consts::FRAC_PI_2).abs() <= 1.0e-12)
            || ((*start_angle_rad + std::f64::consts::FRAC_PI_2).abs() <= 1.0e-12
                && end_angle_rad.abs() <= 1.0e-12)
    )));
    assert!(recovered.max_residual_mm <= 1.0e-12);

    let tampered = split_hemisphere_faces(1.01);
    let edge_faces = edge_face_map(&tampered);
    let context = TopologyContext {
        faces: &tampered,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    assert!(detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5], &tampered, &context).is_none());
    Ok(())
}

fn quarter_fillet_profile() -> Vec<RecoveredProfileCurve> {
    vec![
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [0.0, 1.0],
            end_mm: [1.0, 1.0],
        },
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [1.0, 1.0],
            end_mm: [1.0, 1.4],
        },
        RecoveredProfileCurve::CircleArc {
            source_edge_ids: Vec::new(),
            center_mm: [0.9, 1.4],
            radius_mm: 0.1,
            start_angle_rad: 0.0,
            end_angle_rad: std::f64::consts::FRAC_PI_2,
        },
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [0.9, 1.5],
            end_mm: [0.0, 1.5],
        },
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [0.0, 1.5],
            end_mm: [0.0, 1.0],
        },
    ]
}

#[test]
fn orders_mixed_line_arc_profile_and_rejects_line_arc_crossing() -> anyhow::Result<()> {
    let profile = quarter_fillet_profile();
    let scrambled = vec![
        profile[2].reversed(),
        profile[4].clone(),
        profile[1].clone(),
        profile[3].reversed(),
        profile[0].clone(),
    ];
    let ordered = MeridianProfile::from_recovered(&scrambled, GEOM_TOL_MM)
        .ok_or_else(|| anyhow::anyhow!("expected closed profile"))?
        .into_recovered();
    assert_eq!(ordered.len(), 5);
    assert!(matches!(
        &ordered[2],
        RecoveredProfileCurve::CircleArc {
            center_mm: [0.9, 1.4],
            radius_mm: 0.1,
            start_angle_rad,
            end_angle_rad,
            ..
        } if start_angle_rad.abs() <= 1.0e-12
            && (*end_angle_rad - std::f64::consts::FRAC_PI_2).abs() <= 1.0e-12
    ));

    let crossing = RecoveredProfileCurve::Line {
        source_edge_ids: Vec::new(),
        start_mm: [0.95, 1.2],
        end_mm: [0.95, 1.6],
    };
    assert!(meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::from_recovered(&crossing)
            .ok_or_else(|| anyhow::anyhow!("expected line primitive"))?,
        &MeridianCurve::from_recovered(&profile[2])
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        false,
        GEOM_TOL_MM,
    ));
    Ok(())
}

#[test]
fn orders_two_arc_profile_and_rejects_arc_crossings() -> anyhow::Result<()> {
    let profile = [
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [0.0, 0.0],
            end_mm: [0.9, 0.0],
        },
        RecoveredProfileCurve::CircleArc {
            source_edge_ids: Vec::new(),
            center_mm: [0.9, 0.1],
            radius_mm: 0.1,
            start_angle_rad: -std::f64::consts::FRAC_PI_2,
            end_angle_rad: 0.0,
        },
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [1.0, 0.1],
            end_mm: [1.0, 0.9],
        },
        RecoveredProfileCurve::CircleArc {
            source_edge_ids: Vec::new(),
            center_mm: [0.9, 0.9],
            radius_mm: 0.1,
            start_angle_rad: 0.0,
            end_angle_rad: std::f64::consts::FRAC_PI_2,
        },
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [0.9, 1.0],
            end_mm: [0.0, 1.0],
        },
        RecoveredProfileCurve::Line {
            source_edge_ids: Vec::new(),
            start_mm: [0.0, 1.0],
            end_mm: [0.0, 0.0],
        },
    ];
    let scrambled = vec![
        profile[3].reversed(),
        profile[0].clone(),
        profile[5].clone(),
        profile[2].reversed(),
        profile[4].clone(),
        profile[1].clone(),
    ];
    let ordered = MeridianProfile::from_recovered(&scrambled, GEOM_TOL_MM)
        .ok_or_else(|| anyhow::anyhow!("expected closed profile"))?
        .into_recovered();
    assert_eq!(ordered.len(), 6);
    assert_eq!(
        ordered
            .iter()
            .filter(|curve| matches!(curve, RecoveredProfileCurve::CircleArc { .. }))
            .count(),
        2
    );

    let first = RecoveredProfileCurve::CircleArc {
        source_edge_ids: Vec::new(),
        center_mm: [0.0, 0.0],
        radius_mm: 1.0,
        start_angle_rad: 0.0,
        end_angle_rad: std::f64::consts::PI,
    };
    let crossing = RecoveredProfileCurve::CircleArc {
        source_edge_ids: Vec::new(),
        center_mm: [1.0, 0.0],
        radius_mm: 1.0,
        start_angle_rad: std::f64::consts::FRAC_PI_2,
        end_angle_rad: 3.0 * std::f64::consts::FRAC_PI_2,
    };
    assert!(meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::from_recovered(&first)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        &MeridianCurve::from_recovered(&crossing)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        false,
        GEOM_TOL_MM,
    ));

    let tangent_a = RecoveredProfileCurve::CircleArc {
        source_edge_ids: Vec::new(),
        center_mm: [0.0, 0.0],
        radius_mm: 1.0,
        start_angle_rad: 0.0,
        end_angle_rad: std::f64::consts::FRAC_PI_2,
    };
    let tangent_b = RecoveredProfileCurve::CircleArc {
        source_edge_ids: Vec::new(),
        center_mm: [0.0, 2.0],
        radius_mm: 1.0,
        start_angle_rad: -std::f64::consts::FRAC_PI_2,
        end_angle_rad: 0.0,
    };
    assert!(!meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::from_recovered(&tangent_a)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        &MeridianCurve::from_recovered(&tangent_b)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        true,
        GEOM_TOL_MM,
    ));

    let coincident_adjacent = RecoveredProfileCurve::CircleArc {
        source_edge_ids: Vec::new(),
        center_mm: [0.0, 0.0],
        radius_mm: 1.0,
        start_angle_rad: std::f64::consts::FRAC_PI_2,
        end_angle_rad: std::f64::consts::PI,
    };
    assert!(!meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::from_recovered(&tangent_a)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        &MeridianCurve::from_recovered(&coincident_adjacent)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        true,
        GEOM_TOL_MM,
    ));

    let coincident_overlap = RecoveredProfileCurve::CircleArc {
        source_edge_ids: Vec::new(),
        center_mm: [0.0, 0.0],
        radius_mm: 1.0,
        start_angle_rad: std::f64::consts::FRAC_PI_4,
        end_angle_rad: std::f64::consts::PI,
    };
    assert!(meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::from_recovered(&tangent_a)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        &MeridianCurve::from_recovered(&coincident_overlap)
            .ok_or_else(|| anyhow::anyhow!("expected arc primitive"))?,
        true,
        GEOM_TOL_MM,
    ));
    Ok(())
}

fn split_quarter_torus_faces(
    torus_major_radius_mm: f64,
    torus_minor_radius_mm: f64,
) -> Vec<FaceInfo> {
    let torus = brep::TorusSupport {
        center_mm: [0.0, 0.0, 1.4],
        axis: [0.0, 0.0, 1.0],
        x_direction: [1.0, 0.0, 0.0],
        major_radius_mm: torus_major_radius_mm,
        minor_radius_mm: torus_minor_radius_mm,
        chart_outward: true,
    };
    let cylinder = brep::CylinderSupport {
        axis_origin_mm: [0.0, 0.0, 0.0],
        axis: [0.0, 0.0, 1.0],
        x_direction: [1.0, 0.0, 0.0],
        radius_mm: torus_major_radius_mm + torus_minor_radius_mm,
    };
    let bottom_plane = brep::PlaneSupport {
        origin_mm: [0.0, 0.0, 1.0],
        normal: [0.0, 0.0, 1.0],
        max_residual_mm: 0.0,
    };
    let top_plane = brep::PlaneSupport {
        origin_mm: [0.0, 0.0, 1.4 + torus_minor_radius_mm],
        normal: [0.0, 0.0, 1.0],
        max_residual_mm: 0.0,
    };

    let outer_radius = torus_major_radius_mm + torus_minor_radius_mm;
    let top_z = 1.4 + torus_minor_radius_mm;
    let top_a = test_circle_edge(
        110,
        [torus_major_radius_mm, 0.0, top_z],
        [-torus_major_radius_mm, 0.0, top_z],
        [0.0, 0.0, top_z],
        [0.0, 0.0, 1.0],
        torus_major_radius_mm,
    );
    let top_b = test_circle_edge(
        111,
        [-torus_major_radius_mm, 0.0, top_z],
        [torus_major_radius_mm, 0.0, top_z],
        [0.0, 0.0, top_z],
        [0.0, 0.0, 1.0],
        torus_major_radius_mm,
    );
    let side_a = test_circle_edge(
        112,
        [outer_radius, 0.0, 1.4],
        [-outer_radius, 0.0, 1.4],
        [0.0, 0.0, 1.4],
        [0.0, 0.0, 1.0],
        outer_radius,
    );
    let side_b = test_circle_edge(
        113,
        [-outer_radius, 0.0, 1.4],
        [outer_radius, 0.0, 1.4],
        [0.0, 0.0, 1.4],
        [0.0, 0.0, 1.0],
        outer_radius,
    );
    let bottom_a = test_circle_edge(
        114,
        [outer_radius, 0.0, 1.0],
        [-outer_radius, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        outer_radius,
    );
    let bottom_b = test_circle_edge(
        115,
        [-outer_radius, 0.0, 1.0],
        [outer_radius, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        outer_radius,
    );
    let torus_seam_positive = test_circle_edge(
        116,
        [torus_major_radius_mm, 0.0, top_z],
        [outer_radius, 0.0, 1.4],
        [torus_major_radius_mm, 0.0, 1.4],
        [0.0, 1.0, 0.0],
        torus_minor_radius_mm,
    );
    let torus_seam_negative = test_circle_edge(
        117,
        [-torus_major_radius_mm, 0.0, top_z],
        [-outer_radius, 0.0, 1.4],
        [-torus_major_radius_mm, 0.0, 1.4],
        [0.0, -1.0, 0.0],
        torus_minor_radius_mm,
    );
    let cylinder_seam_positive =
        test_line_edge(118, [outer_radius, 0.0, 1.0], [outer_radius, 0.0, 1.4]);
    let cylinder_seam_negative =
        test_line_edge(119, [-outer_radius, 0.0, 1.0], [-outer_radius, 0.0, 1.4]);

    vec![
        test_face(
            SurfaceSupport::Torus(torus),
            vec![
                top_a.clone(),
                torus_seam_positive.clone(),
                side_a.clone(),
                torus_seam_negative.clone(),
            ],
        ),
        test_face(
            SurfaceSupport::Torus(torus),
            vec![
                top_b.clone(),
                torus_seam_positive,
                side_b.clone(),
                torus_seam_negative,
            ],
        ),
        test_face(
            SurfaceSupport::Cylinder(cylinder),
            vec![
                bottom_a.clone(),
                cylinder_seam_positive.clone(),
                side_a,
                cylinder_seam_negative.clone(),
            ],
        ),
        test_face(
            SurfaceSupport::Cylinder(cylinder),
            vec![
                bottom_b.clone(),
                cylinder_seam_positive,
                side_b,
                cylinder_seam_negative,
            ],
        ),
        test_face(
            SurfaceSupport::Plane(bottom_plane),
            vec![bottom_a, bottom_b],
        ),
        test_face(SurfaceSupport::Plane(top_plane), vec![top_a, top_b]),
    ]
}

#[test]
fn recovers_split_quarter_torus_fillet_topology() -> anyhow::Result<()> {
    let faces = split_quarter_torus_faces(0.9, 0.1);
    let edge_faces = edge_face_map(&faces);
    let entities = Vec::new();
    let index = HashMap::new();
    let context = TopologyContext {
        faces: &faces,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    let recovered = detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5, 6], &faces, &context)
        .ok_or_else(|| anyhow::anyhow!("expected torus-fillet recovery"))?;
    assert_eq!(recovered.profile_curves.len(), 5);
    assert!(recovered.profile_curves.iter().any(|curve| matches!(
        curve,
        RecoveredProfileCurve::CircleArc {
            center_mm: [0.9, 1.4],
            radius_mm: 0.1,
            start_angle_rad,
            end_angle_rad,
            ..
        } if start_angle_rad.abs() <= 1.0e-12
            && (*end_angle_rad - std::f64::consts::FRAC_PI_2).abs() <= 1.0e-12
    )));

    let spindle = split_quarter_torus_faces(0.05, 0.1);
    let edge_faces = edge_face_map(&spindle);
    let context = TopologyContext {
        faces: &spindle,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    let recovered = detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5, 6], &spindle, &context)
        .ok_or_else(|| anyhow::anyhow!("expected spindle-torus recovery"))?;
    assert!(recovered.profile_curves.iter().any(|curve| matches!(
        curve,
        RecoveredProfileCurve::CircleArc {
            center_mm: [0.05, 1.4],
            radius_mm: 0.1,
            start_angle_rad,
            end_angle_rad,
            ..
        } if start_angle_rad.abs() <= 1.0e-12
            && (*end_angle_rad - std::f64::consts::FRAC_PI_2).abs() <= 1.0e-12
    )));

    let quarter = MeridianCurve::Arc(MeridianArc {
        source_edge_ids: Vec::new(),
        center: [0.05, 0.0],
        radius: 0.1,
        start_angle: 0.0,
        end_angle: std::f64::consts::FRAC_PI_2,
    });
    let half = MeridianCurve::Arc(MeridianArc {
        source_edge_ids: Vec::new(),
        center: [0.05, 0.0],
        radius: 0.1,
        start_angle: 0.0,
        end_angle: std::f64::consts::PI,
    });
    assert!(quarter.min_radius(GEOM_TOL_MM) >= 0.0);
    assert!(half.min_radius(GEOM_TOL_MM) < -GEOM_TOL_MM);

    let mut tampered = split_quarter_torus_faces(0.9, 0.1);
    for face in tampered.iter_mut().take(2) {
        let SurfaceSupport::Torus(mut torus) = face.surface else {
            anyhow::bail!("expected torus test face");
        };
        torus.minor_radius_mm = 0.11;
        face.surface = SurfaceSupport::Torus(torus);
    }
    let edge_faces = edge_face_map(&tampered);
    let context = TopologyContext {
        faces: &tampered,
        edge_faces: &edge_faces,
        entities: &entities,
        index: &index,
    };
    assert!(detect_mixed_curved_revolution(99, &[1, 2, 3, 4, 5, 6], &tampered, &context).is_none());
    Ok(())
}

#[test]
fn solid_surface_signature_reports_native_spherical_cap() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_spherical_cap.step");
    let signatures = crate::detect_solid_surface_signatures_bytes(bytes)?;
    assert_eq!(signatures.len(), 1);
    assert_eq!(signatures[0].face_count, 2);
    assert_eq!(signatures[0].support_counts.get("sphere"), Some(&1));
    assert_eq!(signatures[0].support_counts.get("plane"), Some(&1));
    assert!(signatures[0].closed_two_manifold);
    Ok(())
}

#[test]
fn shared_edge_incidence_preserves_twice_used_periodic_seams() -> anyhow::Result<()> {
    let seam = test_line_edge(41, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    let faces = vec![test_face(
        SurfaceSupport::Other { entity_id: 10 },
        vec![seam.clone(), seam],
    )];
    let edges =
        face_edge_uses(&faces).ok_or_else(|| anyhow::anyhow!("twice-used seam was rejected"))?;
    assert_eq!(edges.get(&41), Some(&vec![0, 0]));
    assert_eq!(closed_shell_edge_uses(&faces), Some(edges));
    Ok(())
}

#[test]
fn shared_edge_incidence_distinguishes_open_interfaces_from_closed_shells() -> anyhow::Result<()> {
    let edge = |id| test_line_edge(id, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    let faces = vec![
        test_face(
            SurfaceSupport::Other { entity_id: 11 },
            vec![edge(1), edge(2)],
        ),
        test_face(
            SurfaceSupport::Other { entity_id: 12 },
            vec![edge(1), edge(3)],
        ),
    ];
    let incidence = face_edge_uses(&faces)
        .ok_or_else(|| anyhow::anyhow!("valid open interface not collected"))?;
    assert_eq!(incidence.get(&1), Some(&vec![0, 1]));
    assert_eq!(incidence.get(&2), Some(&vec![0]));
    assert!(selected_faces_connected(2, &incidence));
    assert!(closed_shell_edge_uses(&faces).is_none());

    let nonmanifold = vec![
        faces[0].clone(),
        faces[1].clone(),
        test_face(SurfaceSupport::Other { entity_id: 13 }, vec![edge(1)]),
    ];
    let incident = face_edge_uses(&nonmanifold)
        .ok_or_else(|| anyhow::anyhow!("nonmanifold source incidence missing"))?;
    assert_eq!(incident.get(&1), Some(&vec![0, 1, 2]));
    assert!(!selected_faces_connected(3, &incident));
    assert!(closed_shell_edge_uses(&nonmanifold).is_none());

    // Connectivity must be based on shared two-use edges, not merely the
    // presence of individually admissible one-use rim edges.
    let disconnected =
        HashMap::from([(1, vec![0, 1]), (2, vec![2, 3]), (3, vec![0]), (4, vec![2])]);
    assert!(!selected_faces_connected(4, &disconnected));
    assert!(!shell_faces_connected(4, &disconnected));
    assert!(!selected_faces_connected(
        2,
        &HashMap::from([(1, vec![0, 2])])
    ));

    let no_bound = [FaceInfo {
        surface: SurfaceSupport::Other { entity_id: 20 },
        loops: vec![],
    }];
    assert!(face_edge_uses(&no_bound).is_none());
    assert!(closed_shell_edge_uses(&no_bound).is_none());
    Ok(())
}

#[test]
fn shared_source_face_loader_matches_surface_signature_census() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_spherical_cap.step");
    let exchange = step_p21::parser::parse(std::str::from_utf8(bytes)?)?;
    let entities = &exchange.data[0].entities;
    let index = build_index(entities);
    let solid_id = entities
        .iter()
        .find(|entity| simple_record(entity).is_some_and(|r| r.name == "MANIFOLD_SOLID_BREP"))
        .map(entity_id)
        .ok_or_else(|| anyhow::anyhow!("spherical cap fixture has no solid"))?;
    let face_ids = brep::solid_face_ids(solid_id, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("spherical cap fixture has no faces"))?;
    let faces = source_face_infos(&face_ids, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("source face loading failed"))?;
    let incidence = closed_shell_edge_uses(&faces)
        .ok_or_else(|| anyhow::anyhow!("closed shell incidence failed"))?;
    let census = detect_solid_surface_signatures_with_index(entities, &index);
    let [signature] = census.as_slice() else {
        anyhow::bail!("unexpected number of solid signatures");
    };
    assert_eq!(signature.face_count, faces.len());
    assert_eq!(signature.unique_edge_count, incidence.len());
    assert_eq!(
        signature.edge_use_count,
        incidence.values().map(Vec::len).sum::<usize>()
    );
    assert!(signature.closed_two_manifold);
    Ok(())
}

#[test]
fn shell_connectivity_rejects_disconnected_two_manifolds() {
    let connected = HashMap::from([(1, vec![0, 1]), (2, vec![1, 2]), (3, vec![2, 0])]);
    assert!(shell_faces_connected(3, &connected));

    let disconnected = HashMap::from([
        (1, vec![0, 1]),
        (2, vec![0, 1]),
        (3, vec![2, 3]),
        (4, vec![2, 3]),
    ]);
    assert!(!shell_faces_connected(4, &disconnected));
}

#[test]
fn closes_axis_touching_step_profile() -> anyhow::Result<()> {
    let profile = MeridianProfile::from_lines(
        [
            Segment2 {
                a: [0.0, 0.0],
                b: [2.0, 0.0],
            },
            Segment2 {
                a: [2.0, 0.0],
                b: [2.0, 3.0],
            },
            Segment2 {
                a: [0.0, 3.0],
                b: [2.0, 3.0],
            },
        ]
        .into_iter()
        .map(|segment| MeridianLine {
            start: segment.a,
            end: segment.b,
        }),
        GEOM_TOL_MM,
    )
    .ok_or_else(|| anyhow::anyhow!("expected recovered test geometry"))?
    .into_recovered();
    assert_eq!(profile.len(), 4);
    let points = profile
        .iter()
        .filter_map(RecoveredProfileCurve::start_point)
        .collect::<Vec<_>>();
    assert!(points.contains(&[0.0, 0.0]));
    assert!(points.contains(&[0.0, 3.0]));
    Ok(())
}

#[test]
fn closes_sloped_frustum_profile() -> anyhow::Result<()> {
    let profile = MeridianProfile::from_lines(
        [
            Segment2 {
                a: [0.0, -1.0],
                b: [2.0, -1.0],
            },
            Segment2 {
                a: [2.0, -1.0],
                b: [1.0, 1.0],
            },
            Segment2 {
                a: [0.0, 1.0],
                b: [1.0, 1.0],
            },
        ]
        .into_iter()
        .map(|segment| MeridianLine {
            start: segment.a,
            end: segment.b,
        }),
        GEOM_TOL_MM,
    )
    .ok_or_else(|| anyhow::anyhow!("expected recovered test geometry"))?
    .into_recovered();
    assert_eq!(profile.len(), 4);
    let points = profile
        .iter()
        .filter_map(RecoveredProfileCurve::start_point)
        .collect::<Vec<_>>();
    assert!(points.contains(&[2.0, -1.0]));
    assert!(points.contains(&[1.0, 1.0]));
    Ok(())
}

#[test]
fn rejects_crossing_sloped_profiles() {
    assert!(meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::Line(MeridianLine {
            start: [1.0, 0.0],
            end: [3.0, 2.0]
        }),
        &MeridianCurve::Line(MeridianLine {
            start: [3.0, 0.0],
            end: [1.0, 2.0]
        }),
        false,
        GEOM_TOL_MM,
    ));
    assert!(!meridian::curves_intersect_away_from_shared_endpoint(
        &MeridianCurve::Line(MeridianLine {
            start: [1.0, 0.0],
            end: [2.0, 1.0]
        }),
        &MeridianCurve::Line(MeridianLine {
            start: [3.0, 0.0],
            end: [4.0, 1.0]
        }),
        false,
        GEOM_TOL_MM,
    ));
}

#[test]
fn recovers_native_conical_frustum_fixture() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_conical_frustum.step");
    let recovered = crate::detect_solid_revolutions_bytes(bytes)?;
    assert_eq!(recovered.len(), 1);
    assert_eq!(
        recovered[0]
            .polygon_points()
            .ok_or_else(|| anyhow::anyhow!("expected polygon profile"))?,
        vec![[0.0, 0.0], [2.0, 0.0], [1.0, 2.0], [0.0, 2.0]]
    );
    assert!(recovered[0].max_residual_mm < 1.0e-9);
    assert_eq!(
        recovered[0].source_tolerance_mm,
        REVOLUTION_SOURCE_SUPPORT_TOL_MM
    );

    let declared_uncertainty = String::from_utf8_lossy(bytes).replacen(
        "LENGTH_MEASURE(1.E-07)",
        "LENGTH_MEASURE(1.E-05)",
        1,
    );
    let recovered_with_declared_uncertainty =
        crate::detect_solid_revolutions_bytes(declared_uncertainty.as_bytes())?;
    assert_eq!(recovered_with_declared_uncertainty.len(), 1);
    assert_eq!(
        recovered_with_declared_uncertainty[0].source_tolerance_mm,
        MAX_REVOLUTION_SOURCE_UNCERTAINTY_MM
    );

    let oversized_uncertainty = String::from_utf8_lossy(bytes).replacen(
        "LENGTH_MEASURE(1.E-07)",
        "LENGTH_MEASURE(1.E-04)",
        1,
    );
    let recovered_with_oversized_uncertainty =
        crate::detect_solid_revolutions_bytes(oversized_uncertainty.as_bytes())?;
    assert_eq!(recovered_with_oversized_uncertainty.len(), 1);
    assert_eq!(
        recovered_with_oversized_uncertainty[0].source_tolerance_mm,
        REVOLUTION_SOURCE_SUPPORT_TOL_MM
    );

    let tampered = String::from_utf8_lossy(bytes).replacen(
        "CONICAL_SURFACE('',#32,2.,0.463647609001)",
        "CONICAL_SURFACE('',#32,2.25,0.463647609001)",
        1,
    );
    assert!(crate::detect_solid_revolutions_bytes(tampered.as_bytes())?.is_empty());

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
        step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
    }
    Ok(())
}

#[test]
fn recovers_native_hollow_conical_frustum_fixture() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_hollow_conical_frustum.step");
    let recovered = crate::detect_solid_revolutions_bytes(bytes)?;
    assert_eq!(recovered.len(), 1);
    assert_eq!(
        recovered[0]
            .polygon_points()
            .ok_or_else(|| anyhow::anyhow!("expected polygon profile"))?,
        vec![[0.5, 2.0], [1.0, 0.0], [3.0, 0.0], [2.0, 2.0]]
    );
    assert!(recovered[0].max_residual_mm < 1.0e-9);

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
        step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
    }
    Ok(())
}

#[test]
fn recovers_native_line_surface_of_revolution_fixture() -> anyhow::Result<()> {
    let bytes =
        include_bytes!("../../validation/fixtures/native_line_surface_of_revolution_frustum.step");
    let recovered = crate::detect_solid_revolutions_bytes(bytes)?;
    assert_eq!(recovered.len(), 1);
    assert_eq!(
        recovered[0]
            .polygon_points()
            .ok_or_else(|| anyhow::anyhow!("expected polygon profile"))?,
        vec![[0.0, 0.0], [2.0, 0.0], [1.0, 2.0], [0.0, 2.0]]
    );
    assert!(recovered[0].max_residual_mm < 1.0e-9);

    let skewed = String::from_utf8_lossy(bytes).replacen(
        "CARTESIAN_POINT('',(2.,-4.898587196589E-16,0.))",
        "CARTESIAN_POINT('',(2.,0.25,0.))",
        1,
    );
    assert!(crate::detect_solid_revolutions_bytes(skewed.as_bytes())?.is_empty());

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
    }
    Ok(())
}

#[test]
fn recovers_native_spherical_cap_fixture() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_spherical_cap.step");
    let recovered = crate::detect_solid_revolutions_bytes(bytes)?;
    assert_eq!(recovered.len(), 1);
    assert_eq!(recovered[0].face_ids.len(), 2);
    assert_eq!(recovered[0].profile_curves.len(), 3);

    let plane_axis = [0.0, 0.073];
    let cap_edge = [0.107568582774, 0.073];
    let pole = [0.0, -0.13];
    assert!(profile_has_line(
        &recovered[0].profile_curves,
        plane_axis,
        cap_edge,
    ));
    assert!(profile_has_arc(
        &recovered[0].profile_curves,
        [0.0, 0.0],
        0.13,
        cap_edge,
        pole,
    ));
    assert!(profile_has_line(
        &recovered[0].profile_curves,
        pole,
        plane_axis,
    ));
    assert!(recovered[0].max_residual_mm < 1.0e-9);

    let malformed = String::from_utf8_lossy(bytes).replacen(
        "CIRCLE('',#26,0.107568582774)",
        "CIRCLE('',#26,0.117568582774)",
        1,
    );
    assert!(crate::detect_solid_revolutions_bytes(malformed.as_bytes())?.is_empty());

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
        step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
    }
    Ok(())
}

#[test]
fn recovers_positive_spherical_cap_fixture() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_spherical_cap_positive.step");
    let recovered = crate::detect_solid_revolutions_bytes(bytes)?;
    assert_eq!(recovered.len(), 1);
    assert_eq!(recovered[0].profile_curves.len(), 3);
    let pole = [0.0, 0.13];
    let arc = recovered[0]
        .profile_curves
        .iter()
        .find(|curve| matches!(curve, RecoveredProfileCurve::CircleArc { .. }))
        .ok_or_else(|| anyhow::anyhow!("expected spherical meridian arc"))?;
    let RecoveredProfileCurve::CircleArc {
        center_mm,
        radius_mm,
        ..
    } = arc
    else {
        unreachable!()
    };
    assert!(distance2(*center_mm, [0.0, 0.0]) < 1.0e-12);
    assert!((*radius_mm - 0.13).abs() < 1.0e-12);
    let arc_start = arc
        .start_point()
        .ok_or_else(|| anyhow::anyhow!("arc start missing"))?;
    let arc_end = arc
        .end_point()
        .ok_or_else(|| anyhow::anyhow!("arc end missing"))?;
    assert!(distance2(arc_start, pole) < 1.0e-12 || distance2(arc_end, pole) < 1.0e-12);
    assert!(recovered[0].max_residual_mm < 1.0e-10);

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
        step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
    }
    Ok(())
}

#[test]
fn recovers_native_ring_torus_fixture() -> anyhow::Result<()> {
    let bytes = include_bytes!("../../validation/fixtures/native_torus.step");
    let recovered = crate::detect_solid_revolutions_bytes(bytes)?;
    assert_eq!(recovered.len(), 1);
    assert_eq!(recovered[0].face_ids.len(), 1);
    assert_eq!(recovered[0].profile_curves.len(), 1);
    let RecoveredProfileCurve::CircleArc {
        center_mm,
        radius_mm,
        start_angle_rad,
        end_angle_rad,
        ..
    } = &recovered[0].profile_curves[0]
    else {
        anyhow::bail!("expected full-circle torus meridian");
    };
    assert!((center_mm[0] - 1.2).abs() < 1.0e-12);
    assert!(center_mm[1].abs() < 1.0e-12);
    assert!((*radius_mm - 0.09).abs() < 1.0e-12);
    assert_eq!(*start_angle_rad, 0.0);
    assert_eq!(*end_angle_rad, std::f64::consts::TAU);

    let malformed =
        String::from_utf8_lossy(bytes).replacen("CIRCLE('',#26,1.29)", "CIRCLE('',#26,1.31)", 1);
    assert!(crate::detect_solid_revolutions_bytes(malformed.as_bytes())?.is_empty());

    #[cfg(feature = "cad-kernel-monstertruck")]
    {
        use crate::cad_kernel::CadKernel;
        let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).geometrically_consistent);
        step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
    }
    Ok(())
}

#[test]
fn closes_hollow_step_profile_and_deduplicates_patches() -> anyhow::Result<()> {
    let mut segments = Vec::new();
    for _ in 0..4 {
        push_unique_segment(
            &mut segments,
            Segment2 {
                a: [1.0, 0.0],
                b: [3.0, 0.0],
            },
        );
        push_unique_segment(
            &mut segments,
            Segment2 {
                a: [3.0, 0.0],
                b: [3.0, 2.0],
            },
        );
        push_unique_segment(
            &mut segments,
            Segment2 {
                a: [1.0, 2.0],
                b: [3.0, 2.0],
            },
        );
        push_unique_segment(
            &mut segments,
            Segment2 {
                a: [1.0, 0.0],
                b: [1.0, 2.0],
            },
        );
    }
    assert_eq!(segments.len(), 4);
    let profile = MeridianProfile::from_lines(
        segments.into_iter().map(|segment| MeridianLine {
            start: segment.a,
            end: segment.b,
        }),
        GEOM_TOL_MM,
    )
    .ok_or_else(|| anyhow::anyhow!("expected closed profile"))?
    .into_recovered();
    assert_eq!(profile.len(), 4);
    assert!(
        profile
            .iter()
            .filter_map(RecoveredProfileCurve::start_point)
            .all(|point| point[0] >= 1.0)
    );
    Ok(())
}

#[test]
fn rejects_branching_or_self_intersecting_profiles() {
    assert!(
        MeridianProfile::from_lines(
            [
                Segment2 {
                    a: [1.0, 0.0],
                    b: [3.0, 0.0]
                },
                Segment2 {
                    a: [3.0, 0.0],
                    b: [3.0, 2.0]
                },
                Segment2 {
                    a: [3.0, 2.0],
                    b: [1.0, 2.0]
                },
                Segment2 {
                    a: [1.0, 2.0],
                    b: [1.0, 0.0]
                },
                Segment2 {
                    a: [2.0, 0.0],
                    b: [2.0, 2.0]
                },
            ]
            .into_iter()
            .map(|segment| MeridianLine {
                start: segment.a,
                end: segment.b
            }),
            GEOM_TOL_MM
        )
        .is_none()
    );
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn manufacturer_full_cubic(
    edge_id: u64,
    expected_curve_id: u64,
    entities: &[step_p21::ast::EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
) -> anyhow::Result<
    monstertruck_geometry::prelude::NurbsCurve<monstertruck_geometry::prelude::Vector4>,
> {
    let edge_params =
        crate::step_entities::named_entity_parameters(edge_id, "EDGE_CURVE", entities, index)
            .ok_or_else(|| anyhow::anyhow!("missing manufacturer EDGE_CURVE #{edge_id}"))?;
    let curve_id = edge_params
        .get(3)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow::anyhow!("edge #{edge_id} has no curve reference"))?;
    assert_eq!(curve_id, expected_curve_id);
    // Reuse the conservative trim certificate from the source comparison API:
    // equal endpoint coordinates alone are insufficient if the basis curve
    // revisits an endpoint or contains a disconnected knot discontinuity.
    crate::source_edge_locus::complete_bspline_edge(
        crate::source_edge_locus::StepEdge {
            id: edge_id,
            entities,
        },
        1.0e-8,
    )
}

#[cfg(feature = "cad-kernel-monstertruck")]
#[test]
fn manufacturer_radial_slot_cubic_plane_deviation_bracket() -> anyhow::Result<()> {
    use monstertruck_geometry::analytic_plane_deviation::bound_nurbs_plane_deviation;
    use monstertruck_geometry::prelude::{Point3, Vector3};

    let source = include_bytes!("../../validation/fixtures/production_radial_slot.step");
    let exchange = step_p21::parser::parse(std::str::from_utf8(source)?)?;
    let entities = &exchange.data[0].entities;
    let index = build_index(entities);
    let curve = manufacturer_full_cubic(93481, 60612, entities, &index)?;

    let solid_id = entities
        .iter()
        .find(|entity| simple_record(entity).is_some_and(|r| r.name == "MANIFOLD_SOLID_BREP"))
        .map(entity_id)
        .ok_or_else(|| anyhow::anyhow!("missing source solid"))?;
    let source_tolerance = source_tolerance_by_representation_item(entities, &index)
        .get(&solid_id)
        .copied()
        .ok_or_else(|| anyhow::anyhow!("missing declared source tolerance"))?;
    let faces = brep::solid_face_ids(solid_id, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("missing source faces"))?;
    let mut checked_faces = 0;
    for face_id in faces {
        let loops = brep::face_loops(face_id, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("malformed face loops"))?;
        if !loops
            .iter()
            .flat_map(|loop_| &loop_.edges)
            .any(|edge| edge.edge_id == 93481)
        {
            continue;
        }
        let surface_id = brep::face_surface(face_id, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("missing surface"))?;
        let brep::SurfaceSupport::Plane(plane) =
            brep::surface_support(surface_id, entities, &index)
        else {
            continue;
        };
        let origin = Point3::new(plane.origin_mm[0], plane.origin_mm[1], plane.origin_mm[2]);
        let normal = Vector3::new(plane.normal[0], plane.normal[1], plane.normal[2]);
        let bracket = bound_nurbs_plane_deviation(&curve, origin, normal, 1.0e-10)
            .ok_or_else(|| anyhow::anyhow!("source cubic plane bound failed"))?;
        eprintln!(
            "source edge #93481 on face #{face_id}: source_tol={source_tolerance:.15e} bounds={bracket:?}"
        );
        assert!(bracket.lower <= bracket.upper);
        assert!(
            bracket.exceeds(source_tolerance),
            "source edge should fail its declared uncertainty: {bracket:?}"
        );
        assert!(
            bracket.exceeds(1.54455e-4),
            "source curve also exceeds the previously measured imported OCCT edge tolerance: {bracket:?}"
        );
        assert!(
            bracket.upper - bracket.lower < 2.0e-9,
            "source spline bound is too loose: {bracket:?}"
        );
        checked_faces += 1;
    }
    assert!(
        checked_faces > 0,
        "edge must belong to at least one plane face"
    );
    Ok(())
}

#[cfg(feature = "cad-kernel-monstertruck")]
#[test]
fn manufacturer_radial_slot_cubics_depart_from_their_own_conical_carriers() -> anyhow::Result<()> {
    use monstertruck_geometry::analytic_radial_deviation::bound_nurbs_radial_deviation;
    use monstertruck_geometry::prelude::{Point3, Vector3};

    let source = include_bytes!("../../validation/fixtures/production_radial_slot.step");
    let exchange = step_p21::parser::parse(std::str::from_utf8(source)?)?;
    let [section] = exchange.data.as_slice() else {
        anyhow::bail!("expected one manufacturer source DATA section");
    };
    let entities = &section.entities;
    let index = build_index(entities);
    let source_tolerance = source_tolerance_by_representation_item(entities, &index)
        .get(&15385)
        .copied()
        .ok_or_else(|| anyhow::anyhow!("source representation tolerance missing"))?;
    assert_eq!(source_tolerance, 1.0e-5);

    // STEP face-edge incidence, not a nearest-cone candidate. Both topological
    // edges cover their entire referenced B-spline parameter range, verified
    // by manufacturer_full_cubic before the whole-domain bound is applied.
    for (edge_id, curve_id, cone_face_id, minimum_witness) in [
        (93481, 60612, 58899, 1.5e-5),
        (120654, 122924, 36791, 8.0e-5),
    ] {
        let curve = manufacturer_full_cubic(edge_id, curve_id, entities, &index)?;
        let loops = brep::face_loops(cone_face_id, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("manufacturer conical face has invalid topology"))?;
        assert!(
            loops
                .iter()
                .flat_map(|boundary| &boundary.edges)
                .any(|edge| edge.edge_id == edge_id)
        );
        let surface_id = brep::face_surface(cone_face_id, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("missing authored conical surface"))?;
        let brep::SurfaceSupport::Cone(cone) = brep::surface_support(surface_id, entities, &index)
        else {
            anyhow::bail!("incident source face is not a cone");
        };
        let origin = Point3::new(
            cone.reference_origin_mm[0],
            cone.reference_origin_mm[1],
            cone.reference_origin_mm[2],
        );
        let axis = Vector3::new(cone.axis[0], cone.axis[1], cone.axis[2]);
        let bounds = bound_nurbs_radial_deviation(
            &curve,
            origin,
            axis,
            cone.reference_radius_mm,
            cone.semi_angle_rad.tan(),
            1.0e-9,
        )
        .ok_or_else(|| anyhow::anyhow!("whole-source conical radial bracket unsupported"))?;
        eprintln!(
            "manufacturer conical edge #{edge_id} on face #{cone_face_id}: source_tol={source_tolerance:.13e} bracket={bounds:?}"
        );
        assert!(
            bounds.exceeds(source_tolerance),
            "manufacturer cubic diverges from authored cone above source tolerance: {bounds:?}"
        );
        assert!(
            bounds.lower > minimum_witness,
            "source 3D curve must retain an interior off-cone witness: {bounds:?}"
        );
        assert!(
            bounds.upper - bounds.lower < 1.0e-6,
            "global radial upper bound remains excessively loose: {bounds:?}"
        );
    }
    Ok(())
}
