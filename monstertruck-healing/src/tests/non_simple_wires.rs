use super::*;
use monstertruck_topology::Shell;
use monstertruck_topology::repair::{
    split_non_simple_compressed_wires, split_pinched_compressed_faces,
};

/// RED reproducer for the cell-8 extraction blocker: a face whose boundary
/// wires are each SIMPLE but share one vertex (a PINCHED face, e.g. the
/// pass-through imprint landing a T-junction vertex in two wires) is
/// rejected by `Face::new` via `Wire::disjoint_wires` -- the same
/// `NotSimpleWire` error an intra-wire revisit produces, which misled the
/// cell-8 diagnosis until probed.
#[test]
fn split_pinched_faces_makes_vertex_sharing_boundaries_extractable() {
    let vertices = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 0.0),
        Point2::new(-1.0, -1.0),
    ];
    let line = |a: Point2, b: Point2| TrimmedCurve::new(Line(a, b), (0.0, 1.0));
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: line(vertices[0], vertices[1]),
        },
        CompressedEdge {
            vertices: (1, 2),
            curve: line(vertices[1], vertices[2]),
        },
        CompressedEdge {
            vertices: (2, 0),
            curve: line(vertices[2], vertices[0]),
        },
        CompressedEdge {
            vertices: (0, 3),
            curve: line(vertices[0], vertices[3]),
        },
        CompressedEdge {
            vertices: (3, 4),
            curve: line(vertices[3], vertices[4]),
        },
        CompressedEdge {
            vertices: (4, 0),
            curve: line(vertices[4], vertices[0]),
        },
    ];
    let edge_use = |index: usize| CompressedEdgeIndex {
        index,
        orientation: true,
    };
    let faces = vec![CompressedFace {
        surface: (),
        orientation: true,
        // Two individually SIMPLE closed triangles sharing vertex 0.
        boundaries: vec![
            vec![edge_use(0), edge_use(1), edge_use(2)],
            vec![edge_use(3), edge_use(4), edge_use(5)],
        ],
    }];
    let mut shell = CompressedShell {
        vertices,
        edges,
        faces,
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };
    // RED half: extraction fails on the DISJOINTNESS rule (both wires are
    // simple in isolation).
    assert!(matches!(
        Shell::extract(shell.clone()),
        Err(monstertruck_topology::errors::Error::NotSimpleWire)
    ));

    split_pinched_compressed_faces(&mut shell);

    // GREEN half: the pinch is resolved into one representable face per
    // vertex-sharing wire, same surface, and the shell extracts.
    assert_eq!(shell.faces.len(), 2);
    assert!(shell.faces.iter().all(|face| face.boundaries.len() == 1));
    let extracted = Shell::extract(shell.clone()).expect("pinch split must extract");
    assert_eq!(extracted.len(), 2);
}

/// A single wire revisiting a vertex mid-loop (a bowtie remnant) is split
/// into two closed loops; the resulting PINCHED face (loops share the
/// revisited vertex) is then resolved by the pinch splitter, and the shell
/// extracts as two faces.
#[test]
fn split_non_simple_wires_then_pinch_split_extracts_bowtie() {
    let vertices = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 0.0),
        Point2::new(-1.0, -1.0),
    ];
    let line = |a: Point2, b: Point2| TrimmedCurve::new(Line(a, b), (0.0, 1.0));
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: line(vertices[0], vertices[1]),
        },
        CompressedEdge {
            vertices: (1, 2),
            curve: line(vertices[1], vertices[2]),
        },
        CompressedEdge {
            vertices: (2, 0),
            curve: line(vertices[2], vertices[0]),
        },
        CompressedEdge {
            vertices: (0, 3),
            curve: line(vertices[0], vertices[3]),
        },
        CompressedEdge {
            vertices: (3, 4),
            curve: line(vertices[3], vertices[4]),
        },
        CompressedEdge {
            vertices: (4, 0),
            curve: line(vertices[4], vertices[0]),
        },
    ];
    let edge_use = |index: usize| CompressedEdgeIndex {
        index,
        orientation: true,
    };
    let faces = vec![CompressedFace {
        surface: (),
        orientation: true,
        // ONE wire traversing both triangles through vertex 0 twice.
        boundaries: vec![vec![
            edge_use(0),
            edge_use(1),
            edge_use(2),
            edge_use(3),
            edge_use(4),
            edge_use(5),
        ]],
    }];
    let mut shell = CompressedShell {
        vertices,
        edges,
        faces,
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };
    assert!(matches!(
        Shell::extract(shell.clone()),
        Err(monstertruck_topology::errors::Error::NotSimpleWire)
    ));

    split_non_simple_compressed_wires(&mut shell);
    split_pinched_compressed_faces(&mut shell);

    assert_eq!(shell.faces.len(), 2);
    let extracted = Shell::extract(shell.clone()).expect("bowtie split must extract");
    assert_eq!(extracted.len(), 2);
}

/// A wire retracing one edge both ways (a zero-area SPIKE from the
/// pass-through imprint) is dropped rather than kept as a boundary.
#[test]
fn split_non_simple_wires_drops_pass_through_spikes() {
    let vertices = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(-1.0, 0.0),
        Point2::new(-1.0, -1.0),
    ];
    let line = |a: Point2, b: Point2| TrimmedCurve::new(Line(a, b), (0.0, 1.0));
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: line(vertices[0], vertices[1]),
        },
        CompressedEdge {
            vertices: (0, 2),
            curve: line(vertices[0], vertices[2]),
        },
        CompressedEdge {
            vertices: (2, 3),
            curve: line(vertices[2], vertices[3]),
        },
        CompressedEdge {
            vertices: (3, 0),
            curve: line(vertices[3], vertices[0]),
        },
    ];
    let faces = vec![CompressedFace {
        surface: (),
        orientation: true,
        // Spike out to vertex 1 and straight back, then the real triangle.
        boundaries: vec![vec![
            CompressedEdgeIndex {
                index: 0,
                orientation: true,
            },
            CompressedEdgeIndex {
                index: 0,
                orientation: false,
            },
            CompressedEdgeIndex {
                index: 1,
                orientation: true,
            },
            CompressedEdgeIndex {
                index: 2,
                orientation: true,
            },
            CompressedEdgeIndex {
                index: 3,
                orientation: true,
            },
        ]],
    }];
    let mut shell = CompressedShell {
        vertices,
        edges,
        faces,
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };
    assert!(Shell::extract(shell.clone()).is_err());

    split_non_simple_compressed_wires(&mut shell);

    assert_eq!(shell.faces.len(), 1);
    assert_eq!(shell.faces[0].boundaries.len(), 1);
    assert_eq!(shell.faces[0].boundaries[0].len(), 3);
    assert!(Shell::extract(shell.clone()).is_ok());
}

/// 3-D coincidence alone is not enough evidence to delete a distinct-edge
/// digon. The same topology describes a legitimate periodic seam, so the plain
/// shell path must fail closed and preserve the subloop when no face-local trim
/// geometry exists to distinguish the cases.
#[test]
fn split_non_simple_wires_preserves_geometrically_coincident_digon_without_trim_proof() {
    let vertices = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(-1.0, -1.0, 0.0),
    ];
    let line = |a: Point3, b: Point3| TrimmedCurve::new(Line(a, b), (0.0, 1.0));
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: line(vertices[0], vertices[1]),
        },
        CompressedEdge {
            vertices: (1, 0),
            curve: line(vertices[1], vertices[0]),
        },
        CompressedEdge {
            vertices: (0, 2),
            curve: line(vertices[0], vertices[2]),
        },
        CompressedEdge {
            vertices: (2, 3),
            curve: line(vertices[2], vertices[3]),
        },
        CompressedEdge {
            vertices: (3, 0),
            curve: line(vertices[3], vertices[0]),
        },
    ];
    let edge_use = |index: usize| CompressedEdgeIndex {
        index,
        orientation: true,
    };
    let mut shell = CompressedShell {
        vertices,
        edges,
        faces: vec![CompressedFace {
            surface: (),
            orientation: true,
            boundaries: vec![vec![
                edge_use(0),
                edge_use(1),
                edge_use(2),
                edge_use(3),
                edge_use(4),
            ]],
        }],
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };

    split_non_simple_compressed_wires(&mut shell);
    assert_eq!(shell.faces[0].boundaries.len(), 2);
    let mut lengths = shell.faces[0]
        .boundaries
        .iter()
        .map(Vec::len)
        .collect::<Vec<_>>();
    lengths.sort_unstable();
    assert_eq!(lengths, vec![2, 3]);
    assert!(Shell::extract(shell).is_err());
}

#[test]
fn trimmed_geometric_slit_cleanup_requires_face_local_retrace() {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let edges = vec![
        CompressedEdge {
            vertices: (0, 1),
            curve: Line(p0, p1),
        },
        CompressedEdge {
            vertices: (1, 0),
            curve: Line(p1, p0),
        },
    ];
    let edge_use = |index: usize, trim_curve: Line<Point2>| CompressedEdgeUse {
        index,
        orientation: true,
        trim_curve: Some(trim_curve),
    };
    let make_shell = |left_trim: Line<Point2>, right_trim: Line<Point2>| CompressedTrimmedShell {
        vertices: vec![p0, p1],
        edges: edges.clone(),
        faces: vec![CompressedTrimmedFace {
            surface: (),
            orientation: true,
            boundaries: vec![vec![edge_use(0, left_trim), edge_use(1, right_trim)]],
        }],
    };

    // A true zero-area slit retraces one face-local UV span exactly backwards.
    let mut slit = make_shell(
        Line(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
        Line(Point2::new(1.0, 0.0), Point2::new(0.0, 0.0)),
    );
    remove_geometric_trimmed_slit_boundaries(&mut slit, KERNEL_EPSILON);
    assert!(slit.faces.is_empty());

    // A periodic seam is the same 3-D span but two different UV loci. It must
    // remain topology even though its edge curves are geometrically coincident.
    let mut seam = make_shell(
        Line(Point2::new(0.0, 0.0), Point2::new(0.0, 1.0)),
        Line(
            Point2::new(std::f64::consts::TAU, 1.0),
            Point2::new(std::f64::consts::TAU, 0.0),
        ),
    );
    remove_geometric_trimmed_slit_boundaries(&mut seam, KERNEL_EPSILON);
    assert_eq!(seam.faces.len(), 1);
    assert_eq!(seam.faces[0].boundaries.len(), 1);
    assert_eq!(seam.faces[0].boundaries[0].len(), 2);
}

/// A trimmed periodic-surface face presented as a single boundary wire that
/// walks its two cap circles joined by a doubled SEAM edge -- one edge index
/// used twice, once with each orientation -- revisits the seam-endpoint
/// vertices, so `TrimmedShell::try_from` refuses it with `NotSimpleWire` (the
/// boxy / io1 solidify refusal, spec 007 C3). `split_seam_faces_trimmed`
/// re-forms it as its two vertex-disjoint cap loops (dropping the redundant
/// seam uses, preserving each surviving arc's trim), and it extracts.
#[test]
fn split_seam_faces_resolves_doubled_seam_cylinder() {
    // 0 = seam-bottom, 1 = seam-top, 2 = antipode-bottom, 3 = antipode-top.
    let vertices = vec![(), (), (), ()];
    let edge = |vertices: (usize, usize)| CompressedEdge {
        vertices,
        curve: (),
    };
    let edges = vec![
        edge((0, 2)), // bottom half-arc a
        edge((2, 0)), // bottom half-arc b
        edge((0, 1)), // seam
        edge((1, 3)), // top half-arc a
        edge((3, 1)), // top half-arc b
    ];
    let eu = |index, orientation| CompressedEdgeUse::<()> {
        index,
        orientation,
        trim_curve: None,
    };
    let faces = vec![CompressedTrimmedFace {
        surface: (),
        orientation: true,
        // [bottom_a, bottom_b, SEAM, top_a, top_b, SEAM^-1]
        boundaries: vec![vec![
            eu(0, true),
            eu(1, true),
            eu(2, true),
            eu(3, true),
            eu(4, true),
            eu(2, false),
        ]],
    }];
    let mut shell = CompressedTrimmedShell {
        vertices,
        edges,
        faces,
    };

    // RED: the doubled-seam wire revisits vertices 0 and 1 -> NotSimpleWire.
    assert!(matches!(
        TrimmedShell::try_from(shell.clone()),
        Err(monstertruck_topology::errors::Error::NotSimpleWire)
    ));

    let splits = split_seam_faces_trimmed(&mut shell);
    assert_eq!(splits, 1);

    // GREEN: one face carrying the two vertex-disjoint cap loops; seam dropped.
    assert_eq!(shell.faces.len(), 1);
    assert_eq!(shell.faces[0].boundaries.len(), 2);
    let loop_vertices: Vec<std::collections::BTreeSet<usize>> = shell.faces[0]
        .boundaries
        .iter()
        .map(|wire| {
            wire.iter()
                .flat_map(|edge_use| {
                    let (a, b) = shell.edges[edge_use.index].vertices;
                    [a, b]
                })
                .collect()
        })
        .collect();
    assert!(
        loop_vertices[0].is_disjoint(&loop_vertices[1]),
        "cap loops must be vertex-disjoint",
    );
    assert!(TrimmedShell::try_from(shell).is_ok());
}

#[test]
fn geometric_slit_cleanup_preserves_noncoincident_digon() {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let straight = BsplineCurve::new(KnotVector::bezier_knot(1), vec![p0, p1])
        .expect("valid straight test curve");
    let bowed = BsplineCurve::new(
        KnotVector::bezier_knot(2),
        vec![p1, Point3::new(0.5, 0.25, 0.0), p0],
    )
    .expect("valid bowed test curve");
    let shell = CompressedShell {
        vertices: vec![p0, p1],
        edges: vec![
            CompressedEdge {
                vertices: (0, 1),
                curve: straight,
            },
            CompressedEdge {
                vertices: (1, 0),
                curve: bowed,
            },
        ],
        faces: vec![CompressedFace {
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
                },
            ]],
        }],
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };

    assert!(!geometric_slit(
        &shell.edges,
        &shell.faces[0].boundaries[0],
        1.0e-4,
    ));

    assert_eq!(shell.faces.len(), 1);
    assert_eq!(shell.faces[0].boundaries.len(), 1);
    assert_eq!(shell.faces[0].boundaries[0].len(), 2);
}
