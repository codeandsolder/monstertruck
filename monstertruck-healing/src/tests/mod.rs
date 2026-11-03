use super::*;

mod closed_edges;
mod closed_faces_cylinders;
mod closed_faces_holes;
mod non_simple_wires;
mod shell_orientation;

fn sp<S>(surface: &S, p: Point3, hint: Option<(f64, f64)>) -> Option<(f64, f64)>
where S: SearchParameter<SurfaceParameter, Point = Point3> {
    surface.search_parameter(p, hint, 10)
}

#[test]
#[cfg(feature = "step-test")]
fn step_import() {
    use monstertruck_io::step::load::*;
    const STEP_DIRECTORY: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/step/");
    const STEP_FILES: &[&str] = &[
        "occt-cylinder.step",
        "occt-cone.step",
        "abc-0006.step",
        "abc-0008.step",
        "abc-0035.step",
    ];

    STEP_FILES.iter().for_each(|file_name| {
        println!("{file_name}");
        let path = [STEP_DIRECTORY, file_name].concat();
        let step_string = std::fs::read_to_string(path).unwrap();
        let table = Table::from_step(&step_string).unwrap();
        table.shell.values().for_each(|step_shell| {
            let mut cshell = table.to_compressed_shell(step_shell).unwrap();
            cshell.split_closed_edges_and_faces(0.05);
            monstertruck_topology::Shell::extract(cshell).unwrap();
        });
    });
}

#[test]
fn pass_through_dedup_compacts_edge_stable_ids_with_canonical_edge() {
    use monstertruck_geometry::prelude::Line;
    use monstertruck_topology::{
        StableId,
        compress::{
            CompressedEdge, CompressedEdgeUse, CompressedTrimmedFace, CompressedTrimmedShell,
        },
    };

    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let mut shell = CompressedTrimmedShell {
        vertices: vec![p0, p1],
        edges: vec![
            CompressedEdge {
                vertices: (0, 1),
                curve: Line(p0, p1),
            },
            CompressedEdge {
                vertices: (0, 1),
                curve: Line(p0, p1),
            },
        ],
        faces: vec![
            CompressedTrimmedFace {
                boundaries: vec![vec![CompressedEdgeUse {
                    index: 0,
                    orientation: true,
                    trim_curve: None::<()>,
                }]],
                orientation: true,
                surface: (),
            },
            CompressedTrimmedFace {
                boundaries: vec![vec![CompressedEdgeUse {
                    index: 1,
                    orientation: false,
                    trim_curve: None,
                }]],
                orientation: true,
                surface: (),
            },
        ],
        vertex_stable_ids: Some(vec![StableId::new(1), StableId::new(2)]),
        edge_stable_ids: Some(vec![StableId::new(11), StableId::new(12)]),
        face_stable_ids: Some(vec![StableId::new(21), StableId::new(22)]),
    };

    dedup_pass_through_trimmed(&mut shell, 1.0e-9);

    assert_eq!(shell.edges.len(), 1);
    assert_eq!(shell.edge_stable_ids, Some(vec![StableId::new(11)]));
    assert!(
        shell
            .faces
            .iter()
            .flat_map(|face| face.boundaries.iter())
            .flatten()
            .all(|edge_use| edge_use.index == 0)
    );
    assert_eq!(shell.vertex_stable_ids.as_ref().map(Vec::len), Some(2));
    assert_eq!(shell.face_stable_ids.as_ref().map(Vec::len), Some(2));
}
