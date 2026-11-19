//! Transactional fillet shell editing: commit on success, preserve on failure.
use super::*;
use crate::ops::edit_shell_transactionally;

#[test]
fn partial_topology_mutations_are_not_visible_after_failure() {
    let (mut shell, _, _) = build_box_shell();
    let original_faces = shell.face_iter().map(Face::id).collect::<Vec<_>>();
    let original_orientations = shell.face_iter().map(Face::orientation).collect::<Vec<_>>();
    let original_edges = shell.edge_iter().map(|edge| edge.id()).collect::<Vec<_>>();

    let error = edit_shell_transactionally(&mut shell, |candidate| {
        candidate[0].invert();
        candidate.remove(1);
        candidate.pop();
        Err(FilletError::GeometryFailed {
            context: "intentional late topology construction failure",
        })
    });
    assert!(matches!(error, Err(FilletError::GeometryFailed { .. })));
    assert_eq!(
        shell.face_iter().map(Face::id).collect::<Vec<_>>(),
        original_faces
    );
    assert_eq!(
        shell.face_iter().map(Face::orientation).collect::<Vec<_>>(),
        original_orientations,
    );
    assert_eq!(
        shell.edge_iter().map(|edge| edge.id()).collect::<Vec<_>>(),
        original_edges
    );
}

#[test]
fn successful_topology_transaction_commits_once() {
    let (mut shell, _, _) = build_box_shell();
    let initial_faces = shell.len();
    let initial_first_face_orientation = shell[0].orientation();
    let first_face_id = shell[0].id();

    edit_shell_transactionally(&mut shell, |candidate| {
        candidate[0].invert();
        candidate.pop();
        Ok(())
    })
    .unwrap();

    assert_eq!(shell.len(), initial_faces - 1);
    assert_eq!(shell[0].orientation(), !initial_first_face_orientation);
    assert_eq!(shell[0].id(), first_face_id);
}

fn late_chain_failure_fixture() -> (Shell, [EdgeId; 5], FilletOptions) {
    // The three connected top edges form a longer open chain, so it is
    // processed before the separate two-edge closed chain below.
    let (mut shell, box_edges, _) = build_5face_box();
    let radius = FilletOptions {
        radius: RadiusSpec::Variable(Box::new(|t| 0.15 + 0.05 * t)),
        ..Default::default()
    };
    let open_edges = [box_edges[0].id(), box_edges[1].id(), box_edges[2].id()];
    // A second component has two distinct manifold edges forming a closed
    // wire. The same variable radius is invalid there (f(0) != f(1)).
    let front = Vertex::new(Point3::new(10.0, 0.0, 0.0));
    let back = Vertex::new(Point3::new(11.0, 0.0, 0.0));
    let make_edge = |v0: &Vertex, v1: &Vertex| {
        Edge::new_unchecked(
            v0,
            v1,
            NurbsCurve::from(BsplineCurve::new_unchecked(
                KnotVector::bezier_knot(1),
                vec![v0.point(), v1.point()],
            ))
            .into(),
        )
    };
    let forward = make_edge(&front, &back);
    let backward = make_edge(&back, &front);
    let knots = KnotVector::bezier_knot(1);
    let surface = BsplineSurface::new_unchecked(
        (knots.clone(), knots),
        vec![
            vec![front.point(), front.point() + Vector3::new(0.0, 1.0, 0.0)],
            vec![back.point(), back.point() + Vector3::new(0.0, 1.0, 0.0)],
        ],
    );
    let closing_face = Face::new_unchecked(
        vec![vec![forward.clone(), backward.clone()].into()],
        surface.into(),
    );
    shell.push(closing_face.clone());
    shell.push(closing_face.inverse());

    let selected_edges = [
        open_edges[0],
        open_edges[1],
        open_edges[2],
        forward.id(),
        backward.id(),
    ];
    (shell, selected_edges, radius)
}

fn assert_shell_topology_unchanged(actual: &Shell, original: &Shell) {
    assert_eq!(
        actual.face_iter().map(Face::id).collect::<Vec<_>>(),
        original.face_iter().map(Face::id).collect::<Vec<_>>(),
    );
    assert_eq!(
        actual
            .face_iter()
            .map(Face::orientation)
            .collect::<Vec<_>>(),
        original
            .face_iter()
            .map(Face::orientation)
            .collect::<Vec<_>>(),
    );
    assert_eq!(
        actual.edge_iter().map(|edge| edge.id()).collect::<Vec<_>>(),
        original
            .edge_iter()
            .map(|edge| edge.id())
            .collect::<Vec<_>>(),
    );
}

#[test]
fn fillet_request_rolls_back_earlier_chains_after_late_hard_failure() {
    let (mut shell, selected_edges, radius) = late_chain_failure_fixture();
    let original = shell.clone();
    let mut first_chain_only = shell.clone();
    fillet_edges(&mut first_chain_only, &selected_edges[..3], Some(&radius)).unwrap();
    assert!(
        first_chain_only.len() > shell.len(),
        "the initial chain must actually change the shell",
    );
    let result = fillet_edges(&mut shell, &selected_edges, Some(&radius));
    assert!(
        matches!(result, Err(FilletError::VariableRadiusUnsupported)),
        "{result:?}"
    );
    assert_shell_topology_unchanged(&shell, &original);
}

#[test]
fn generic_fillet_request_preserves_external_shell_after_late_hard_failure() {
    let (mut shell, selected_ids, radius) = late_chain_failure_fixture();
    let original = shell.clone();
    let selected_edges = selected_ids
        .into_iter()
        .map(|id| shell.edge_iter().find(|edge| edge.id() == id).unwrap())
        .collect::<Vec<_>>();
    let result = fillet_edges_generic(&mut shell, &selected_edges, Some(&radius));
    assert!(
        matches!(result, Err(FilletError::VariableRadiusUnsupported)),
        "{result:?}"
    );
    assert_shell_topology_unchanged(&shell, &original);
}
