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
