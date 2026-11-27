//! External-crate API regression: a native proof snapshot is read-only,
//! retains its original numerical thresholds, and exposes its entire graph.

use monstertruck_modeling::{Point3, Vector3, builder, certify_native_analytic_boundary_incidence};

#[test]
fn downstream_reader_can_inspect_complete_checked_boundary_graph() {
    let vertex = builder::vertex(Point3::new(-0.5, -0.5, -0.5));
    let edge = builder::extrude(&vertex, Vector3::unit_x());
    let face = builder::extrude(&edge, Vector3::unit_y());
    let solid = builder::extrude(&face, Vector3::unit_z());

    let linear = 1.0e-7;
    let angular = 1.0e-9;
    let result = certify_native_analytic_boundary_incidence(&solid, linear, angular);
    assert!(result.is_ok());
    let Ok(proof) = result else { return };
    assert_eq!(proof.linear_tolerance(), linear);
    assert_eq!(proof.angular_tolerance(), angular);

    let shells = proof.shells();
    assert_eq!(shells.len(), 1);
    let shell = &shells[0];
    assert_eq!(shell.faces().len(), 6);
    assert_eq!(shell.vertices().len(), 8);
    assert_eq!(shell.face_wires().len(), 6);
    assert_eq!(shell.edges().len(), 12);

    for (face_index, wires) in shell.face_wires().iter().enumerate() {
        let face = shell.faces()[face_index];
        let _carrier = face.carrier();
        let _normal_sense = face.face_orientation();
        for (wire_index, wire) in wires.iter().enumerate() {
            assert!(!wire.edges().is_empty());
            for (edge_index, use_) in wire.edges().iter().copied().enumerate() {
                let next = wire.edges()[(edge_index + 1) % wire.edges().len()];
                assert_eq!(use_.face_index(), face_index);
                assert_eq!(use_.wire_index(), wire_index);
                assert_eq!(use_.wire_edge_index(), edge_index);
                assert_eq!(use_.end_vertex(), next.start_vertex());
                let start = shell.vertices()[use_.start_vertex()].point();
                let end = shell.vertices()[use_.end_vertex()].point();
                for coordinate in [start.x, start.y, start.z, end.x, end.y, end.z] {
                    assert!(coordinate.is_finite());
                }
            }
        }
    }
    for edge in shell.edges().iter().copied() {
        let [first, second] = edge.uses();
        assert_eq!(first.start_vertex(), second.end_vertex());
        assert_eq!(first.end_vertex(), second.start_vertex());
        assert!(!edge.is_seam());
        assert!(!edge.same_oriented_support());
        assert!(edge.finite_segment_locus());
    }
}
