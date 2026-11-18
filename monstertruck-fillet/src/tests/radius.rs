use monstertruck_meshing::prelude::*;

// `FilletError` is deliberately NOT imported here: the tests below name it as
// `super::FilletError`, and importing it would make that qualification
// redundant.
use super::{
    FilletOptions, RadiusSpec, build_5face_box, build_box_shell, build_stretched_open_box,
    fillet_along_wire, fillet_closed_top_loop, fillet_edges, prepare_stretched_open_wire,
};

/// Variable-radius fillet along a closed wire (radius varies 0.15..0.20, f(0) ≈ f(1)).
#[test]
fn variable_radius_closed_wire() {
    // Closed-wire continuity requires equal endpoint radii.
    fillet_closed_top_loop(&FilletOptions {
        radius: RadiusSpec::Variable(Box::new(|t| 0.15 + 0.05 * (std::f64::consts::PI * t).sin())),
        ..Default::default()
    });
}

// ---------------------------------------------------------------------------
// Phase 6c: Variable radius on open wires
// ---------------------------------------------------------------------------

/// Variable radius on an open wire should succeed (no f(0)≈f(1) constraint).
#[test]
fn variable_radius_open_wire() {
    // Reuse fillet_semi_cube topology: 4-face open box, 2-edge open wire.
    let (mut shell, edge, v) = build_stretched_open_box();

    let opts = FilletOptions {
        radius: RadiusSpec::Constant(0.4),
        ..Default::default()
    };
    let boundary = prepare_stretched_open_wire(&mut shell, &edge, &opts);
    assert_eq!(boundary.front_vertex().unwrap(), &v[0]);
    assert!(!boundary.is_closed());

    // Variable radius where f(0)=0.1, f(1)=0.3 -- NOT equal, would fail on closed wire.
    let var_opts = FilletOptions {
        radius: RadiusSpec::Variable(Box::new(|t| 0.1 + 0.2 * t)),
        ..Default::default()
    };
    fillet_along_wire(&mut shell, &boundary, &var_opts).unwrap();

    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

// ---------------------------------------------------------------------------
// Phase 6b: Per-edge radius
// ---------------------------------------------------------------------------

/// Per-edge radius with two edges having different radii.
#[test]
fn per_edge_radius_two_edges() {
    let (mut shell, edge, _) = build_5face_box();

    let initial_face_count = shell.len();

    // Two independent edges with different radii.
    let params = FilletOptions {
        radius: RadiusSpec::PerEdge(vec![0.3, 0.15]),
        ..Default::default()
    };
    fillet_edges(&mut shell, &[edge[5].id(), edge[7].id()], Some(&params)).unwrap();

    assert!(
        shell.len() >= initial_face_count + 2,
        "expected at least 2 new fillet faces, got {} total (was {})",
        shell.len(),
        initial_face_count
    );
    let _poly = shell.robust_triangulation(0.001).to_polygon();
}

/// Per-edge radius count mismatch → PerEdgeRadiusMismatch error.
#[test]
fn per_edge_radius_mismatch() {
    let (mut shell, edge, _) = build_box_shell();

    // Provide 1 radius for 2 edges → mismatch.
    let params = FilletOptions {
        radius: RadiusSpec::PerEdge(vec![0.3]),
        ..Default::default()
    };
    let result = fillet_edges(&mut shell, &[edge[5].id(), edge[6].id()], Some(&params));
    assert!(
        matches!(
            result,
            Err(super::FilletError::PerEdgeRadiusMismatch {
                given: 1,
                expected: 2
            })
        ),
        "expected PerEdgeRadiusMismatch, got: {result:?}"
    );
}

/// Per-edge radius where one edge is too short → DegenerateEdge.
#[test]
fn per_edge_radius_degenerate() {
    let (mut shell, edge, _) = build_box_shell();

    // edge[5] length ~1.0, radius 0.15 → ok (2*0.15=0.3 < 1.0).
    // edge[6] length ~1.0, radius 0.6 → too big (2*0.6=1.2 > 1.0).
    let params = FilletOptions {
        radius: RadiusSpec::PerEdge(vec![0.15, 0.6]),
        ..Default::default()
    };
    let result = fillet_edges(&mut shell, &[edge[5].id(), edge[6].id()], Some(&params));
    assert!(
        matches!(result, Err(super::FilletError::DegenerateEdge)),
        "expected DegenerateEdge, got: {result:?}"
    );
}
