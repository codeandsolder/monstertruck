use super::*;
use monstertruck_topology::Shell;
use std::f64::consts::PI;

#[test]
fn test_split_closed_edges() {
    let vertices = vec![Point2::new(1.0, 0.0)];
    let curve = TrimmedCurve::new(UnitCircle::<Point2>::new(), (0.0, 2.0 * PI));
    let edges = vec![CompressedEdge {
        vertices: (0, 0),
        curve,
    }];
    let faces = vec![
        CompressedFace {
            surface: (),
            orientation: true,
            boundaries: vec![vec![CompressedEdgeIndex {
                index: 0,
                orientation: true,
            }]],
        },
        CompressedFace {
            surface: (),
            orientation: false,
            boundaries: vec![vec![CompressedEdgeIndex {
                index: 0,
                orientation: false,
            }]],
        },
    ];
    let mut shell = CompressedShell {
        vertices,
        edges,
        faces,
        vertex_stable_ids: Some(vec![StableId::new(11)]),
        edge_stable_ids: Some(vec![StableId::new(21)]),
        face_stable_ids: Some(vec![StableId::new(31), StableId::new(32)]),
    };
    assert!(Shell::extract(shell.clone()).is_err());

    split_closed_edges(&mut shell);
    assert!(Shell::extract(shell.clone()).is_ok());

    let CompressedShell {
        vertices,
        edges,
        faces,
        ..
    } = &shell;
    assert_eq!(vertices.len(), 2);
    assert_near!(vertices[0], Point2::new(1.0, 0.0));
    assert_near!(vertices[1], Point2::new(-1.0, 0.0));

    assert_eq!(edges.len(), 2);
    assert_eq!(edges[0].vertices, (0, 1));
    assert_eq!(edges[1].vertices, (1, 0));
    assert_eq!(
        shell.vertex_stable_ids.as_deref(),
        Some([StableId::new(11), StableId::UNASSIGNED].as_slice())
    );
    assert_eq!(
        shell.edge_stable_ids.as_deref(),
        Some([StableId::new(21), StableId::UNASSIGNED].as_slice())
    );
    assert_eq!(
        shell.face_stable_ids.as_deref(),
        Some([StableId::new(31), StableId::new(32)].as_slice())
    );

    assert_eq!(
        *faces,
        vec![
            CompressedFace {
                surface: (),
                orientation: true,
                boundaries: vec![vec![
                    CompressedEdgeIndex {
                        index: 0,
                        orientation: true,
                    },
                    CompressedEdgeIndex {
                        index: 1,
                        orientation: true,
                    }
                ]],
            },
            CompressedFace {
                surface: (),
                orientation: false,
                boundaries: vec![vec![
                    CompressedEdgeIndex {
                        index: 1,
                        orientation: false,
                    },
                    CompressedEdgeIndex {
                        index: 0,
                        orientation: false,
                    }
                ]],
            },
        ]
    );
}
