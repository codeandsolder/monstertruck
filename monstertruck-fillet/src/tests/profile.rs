use super::*;

/// Exercise the same two-edge open-wire geometry for the supplied profile,
/// retaining the 0.4 carrier-filleting and 0.2 final-wire radii.
fn exercise_open_wire_profile(profile: FilletProfile) {
    let (mut shell, edges, _) = build_stretched_open_box();
    let initial = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        profile: profile.clone(),
        ..Default::default()
    };
    let boundary = prepare_stretched_open_wire(&mut shell, &edges, &initial);
    let final_options = FilletOptions {
        radius: RadiusSpec::Constant(0.2),
        profile,
        ..Default::default()
    };
    fillet_along_wire(&mut shell, &boundary, &final_options).unwrap();
    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

// ---------------------------------------------------------------------------
// Chamfer tests
// ---------------------------------------------------------------------------

/// Single-edge chamfer on a 2-face shell.
#[test]
fn chamfer_single_edge() {
    let (mut shell, edge, _) = build_box_shell();
    let initial_face_count = shell.len();

    let params = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        profile: FilletProfile::Chamfer,
        ..Default::default()
    };
    fillet_edges(&mut shell, &[edge[5].id()], Some(&params)).unwrap();

    assert!(shell.len() > initial_face_count);
    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

/// Chamfer along an open wire (same topology as fillet_semi_cube).
#[test]
fn chamfer_semi_cube() { exercise_open_wire_profile(FilletProfile::Chamfer); }

#[test]
fn chamfer_closed_wire() {
    fillet_closed_top_loop(&FilletOptions {
        radius: RadiusSpec::Constant(0.2),
        profile: FilletProfile::Chamfer,
        ..Default::default()
    });
}

// ---------------------------------------------------------------------------
// Ridge tests
// ---------------------------------------------------------------------------

/// Single-edge ridge on a 2-face shell.
#[test]
fn ridge_single_edge() {
    let (mut shell, edge, _) = build_box_shell();
    let initial_face_count = shell.len();

    let params = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        profile: FilletProfile::Ridge,
        ..Default::default()
    };
    fillet_edges(&mut shell, &[edge[5].id()], Some(&params)).unwrap();

    assert!(shell.len() > initial_face_count);
    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

/// Ridge along an open wire (same topology as chamfer_semi_cube).
#[test]
fn ridge_semi_cube() { exercise_open_wire_profile(FilletProfile::Ridge); }

#[test]
fn ridge_closed_wire() {
    fillet_closed_top_loop(&FilletOptions {
        radius: RadiusSpec::Constant(0.2),
        profile: FilletProfile::Ridge,
        ..Default::default()
    });
}

// ---------------------------------------------------------------------------
// Custom profile tests
// ---------------------------------------------------------------------------

/// Custom with linear profile (0,0)→(1,0) -- should behave like chamfer.
#[test]
fn custom_profile_linear() {
    let (mut shell, edge, _) = build_box_shell();
    let initial_face_count = shell.len();

    let profile = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
    );
    let params = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        profile: FilletProfile::Custom(Box::new(profile)),
        ..Default::default()
    };
    fillet_edges(&mut shell, &[edge[5].id()], Some(&params)).unwrap();

    assert!(shell.len() > initial_face_count);
    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

/// Custom with degree-2 bump (0,0)→(0.5,1.0)→(1,0) -- non-trivial shape.
#[test]
fn custom_profile_bump() {
    let (mut shell, edge, _) = build_box_shell();
    let initial_face_count = shell.len();

    let profile = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(2),
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 1.0),
            Point2::new(1.0, 0.0),
        ],
    );
    let params = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        profile: FilletProfile::Custom(Box::new(profile)),
        ..Default::default()
    };
    fillet_edges(&mut shell, &[edge[5].id()], Some(&params)).unwrap();

    assert!(shell.len() > initial_face_count);
    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

#[test]
fn asymmetric_chamfer_uses_independent_carrier_setbacks() {
    let (shell, edge, _) = build_box_shell();
    let params = FilletOptions {
        radius: RadiusSpec::Constant(0.1),
        profile: FilletProfile::AsymmetricChamfer {
            first_setback_mm: 0.2,
            second_setback_mm: 0.35,
        },
        ..Default::default()
    };
    let (_, _, chamfer) = fillet(&shell[0], &shell[1], edge[0].id(), &params).unwrap();
    let boundary = chamfer.absolute_boundaries();
    let first_contact = &boundary[0][0];
    let second_contact = &boundary[0][2];
    for vertex in [first_contact.front(), first_contact.back()] {
        let point = vertex.point();
        assert!((point.y - 0.2).abs() < 1.0e-8);
        assert!((point.z - 1.0).abs() < 1.0e-8);
    }
    for vertex in [second_contact.front(), second_contact.back()] {
        let point = vertex.point();
        assert!(point.y.abs() < 1.0e-8);
        assert!((point.z - 0.65).abs() < 1.0e-8);
    }
}

#[test]
fn asymmetric_chamfer_does_not_use_round_radius_length_guard() {
    let (mut shell, edge, _) = build_box_shell();
    let initial_face_count = shell.len();
    let params = FilletOptions {
        radius: RadiusSpec::Constant(10.0),
        profile: FilletProfile::AsymmetricChamfer {
            first_setback_mm: 0.2,
            second_setback_mm: 0.35,
        },
        ..Default::default()
    };
    fillet_edges(&mut shell, &[edge[0].id()], Some(&params)).unwrap();
    assert!(shell.len() > initial_face_count);
}

#[test]
fn asymmetric_chamfer_rejects_invalid_setback_and_multi_edge_selection() {
    let (mut shell, edge, _) = build_box_shell();
    let invalid = FilletOptions {
        radius: RadiusSpec::Constant(0.1),
        profile: FilletProfile::AsymmetricChamfer {
            first_setback_mm: 0.0,
            second_setback_mm: 0.2,
        },
        ..Default::default()
    };
    assert!(matches!(
        fillet(&shell[0], &shell[1], edge[0].id(), &invalid),
        Err(FilletError::InvalidAsymmetricChamferSetback)
    ));

    let multi = FilletOptions {
        radius: RadiusSpec::Constant(0.1),
        profile: FilletProfile::AsymmetricChamfer {
            first_setback_mm: 0.2,
            second_setback_mm: 0.3,
        },
        ..Default::default()
    };
    assert!(matches!(
        fillet_edges(&mut shell, &[edge[0].id(), edge[2].id()], Some(&multi)),
        Err(FilletError::AsymmetricChamferWireUnsupported)
    ));
}
