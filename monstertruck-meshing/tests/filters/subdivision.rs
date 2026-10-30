use monstertruck_mesh::*;
use monstertruck_meshing::filters::*;

#[test]
fn loop_subdivision() {
    let positions = vec![
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, -1.0),
        Point3::new(1.0, 0.0, 0.0),
    ];
    let tri_faces = vec![
        [0.into(), 1.into(), 2.into()],
        [3.into(), 2.into(), 1.into()],
    ];
    let mut polygon = PolygonMesh::new(
        StandardAttributes {
            positions,
            uv_coords: Vec::new(),
            normals: Vec::new(),
        },
        Faces::from_tri_and_quad_faces(tri_faces, Vec::new()),
    )
    .expect("fixture mesh is valid");
    polygon.loop_subdivision();
    assert_eq!(polygon.positions().len(), 9);
    assert_eq!(
        polygon.positions(),
        &[
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 6.0, 6.0) / 8.0,
            Point3::new(0.0, 6.0, -6.0) / 8.0,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 0.75, 0.0),
            Point3::new(-0.5, 0.5, -0.5),
            Point3::new(-0.5, 0.5, 0.5),
            Point3::new(0.5, 0.5, 0.5),
            Point3::new(0.5, 0.5, -0.5),
        ],
    );
}

#[test]
fn loop_subdivision_leaves_invalid_index_mesh_unchanged() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let faces = Faces::from_tri_and_quad_faces(vec![[0.into(), 1.into(), 3.into()]], Vec::new());
    let mut polygon = PolygonMesh::new_unchecked(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        faces,
    );
    let before_positions = polygon.positions().to_vec();
    let before_faces = polygon.faces().clone();

    polygon.loop_subdivision();

    assert_eq!(polygon.positions().as_slice(), before_positions.as_slice());
    assert_eq!(polygon.faces(), &before_faces);
}
