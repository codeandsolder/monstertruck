use super::*;

type CompressedSolid = monstertruck_topology::compress::CompressedSolid<Point3, Curve, Surface>;

macro_rules! dir ( () => { concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/shape/") });

const SHAPE_JSONS: &[&str] = &[
    concat!(dir!(), "bottle.json"),
    concat!(dir!(), "punched-cube.json"),
    concat!(dir!(), "torus-punched-cube.json"),
    concat!(dir!(), "sphere.json"),
    concat!(dir!(), "torus.json"),
];

fn read_jsons() -> Vec<Vec<u8>> {
    let closure = |json| std::fs::read(json).unwrap();
    SHAPE_JSONS.iter().map(closure).collect()
}

#[test]
fn solid_is_closed() {
    for (i, json) in read_jsons().into_iter().enumerate() {
        let solid: Solid = serde_json::from_reader(json.as_slice()).unwrap();
        let mut poly = solid.triangulation(0.01).to_polygon();
        poly.put_together_same_attrs(TOLERANCE * 2.0)
            .remove_degenerate_faces()
            .remove_unused_attrs();
        assert_eq!(
            poly.shell_condition(),
            ShellCondition::Closed,
            "not closed: file no. {i}"
        );
    }
}

#[test]
fn csolid_is_closed() {
    for (i, json) in read_jsons().into_iter().enumerate() {
        let solid: CompressedSolid = serde_json::from_reader(json.as_slice()).unwrap();
        let mut poly = solid.triangulation(0.01).to_polygon();
        poly.put_together_same_attrs(TOLERANCE * 2.0)
            .remove_degenerate_faces()
            .remove_unused_attrs();
        assert_eq!(
            poly.shell_condition(),
            ShellCondition::Closed,
            "not closed: file no. {i}"
        );
    }
}

fn assert_torus_punched_cube_mesh(mut mesh: PolygonMesh) {
    mesh.put_together_same_attrs(TOLERANCE * 2.0)
        .remove_degenerate_faces()
        .remove_unused_attrs();

    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);

    // The fixture is constructed in `monstertruck-modeling/examples/torus-punched-cube.rs`:
    // a unit cube with a quarter of a torus tube removed.  The torus has
    // major radius 1/2 and minor radius 1/4, so the removed quarter-volume is
    // (1/4) * 2*pi^2*R*r^2 = pi^2/64.  This is an independent geometric
    // oracle; no mesh from OCCT/FreeCAD participates in the expected answer.
    let expected_volume = 1.0 - std::f64::consts::PI.powi(2) / 64.0;
    let volume = mesh.volume();
    assert!(
        (volume - expected_volume).abs() < 3.0e-3,
        "torus-punched cube mesh volume {volume:.9}, expected about {expected_volume:.9}"
    );

    // Every tessellation vertex must lie on one of the six cube planes or on
    // the exact torus surface (sqrt(x^2+y^2)-R)^2 + (z-z0)^2 = r^2.
    // This catches a mesh that has roughly the right volume but wanders off the
    // modeled boundary.
    for &point in mesh.positions() {
        let on_cube = [point.x, point.y, point.z]
            .into_iter()
            .any(|coordinate| coordinate.abs() < 1.0e-7 || (coordinate - 1.0).abs() < 1.0e-7);
        let radial = (point.x * point.x + point.y * point.y).sqrt();
        let torus_distance =
            (((radial - 0.5).powi(2) + (point.z - 0.5).powi(2)).sqrt() - 0.25).abs();
        assert!(
            on_cube || torus_distance < 1.0e-7,
            "mesh vertex {point:?} lies on neither the cube nor the analytic torus"
        );
    }
}

#[test]
fn torus_punched_cube_mesh_matches_analytic_geometry() {
    let jsons = read_jsons();
    let solid: Solid = serde_json::from_slice(jsons[2].as_slice()).unwrap();
    assert_torus_punched_cube_mesh(solid.triangulation(0.005).to_polygon());
}

#[test]
fn compressed_torus_punched_cube_mesh_matches_analytic_geometry() {
    let jsons = read_jsons();
    let solid: CompressedSolid = serde_json::from_slice(jsons[2].as_slice()).unwrap();
    assert_torus_punched_cube_mesh(solid.triangulation(0.005).to_polygon());
}

#[test]
fn large_number_meshing() {
    let json = std::fs::read(concat!(dir!(), "large-torus.json")).unwrap();
    let torus: Solid = serde_json::from_slice(json.as_slice()).unwrap();
    let _ = torus.triangulation(1.0).to_polygon();
}

fn special_cylinder_model() -> Shell {
    let v0 = builder::vertex(Point3::new(0.0, 1.0, 0.0));
    let v1 = builder::vertex(Point3::new(0.0, 1.0, 1.0));
    let v2 = builder::vertex(Point3::new(0.0, -1.0, 0.0));
    let v3 = builder::vertex(Point3::new(0.0, -1.0, 1.0));

    let edge0 = builder::line(&v0, &v1);
    let edge1 = builder::line(&v2, &v3);
    let edge2 = builder::circle_arc(&v0, &v2, Point3::new(-1.0, 0.0, 0.0));
    let edge3 = builder::circle_arc(&v2, &v0, Point3::new(1.0, 0.0, 0.0));
    let edge4 = builder::circle_arc(&v1, &v3, Point3::new(-1.0, 0.0, 1.0));
    let edge5 = builder::circle_arc(&v3, &v1, Point3::new(1.0, 0.0, 1.0));

    let face0 =
        builder::try_attach_plane(&[vec![edge2.inverse(), edge3.inverse()].into()]).unwrap();
    let face1 = builder::try_attach_plane(&[vec![edge4.clone(), edge5.clone()].into()]).unwrap();

    let surface_row = RevolutionSurface::<Curve>::by_revolution(
        Line(Point3::new(1.0, 0.0, 1.0), Point3::new(1.0, 0.0, 0.0)).into(),
        Point3::origin(),
        Vector3::unit_z(),
    );
    let surface: Surface = Processor::new(surface_row).into();

    let face2 = Face::try_new(
        vec![vec![edge2, edge1.clone(), edge4.inverse(), edge0.inverse()].into()],
        surface.clone(),
    )
    .unwrap();
    let face3 = Face::try_new(
        vec![vec![edge3, edge0, edge5.inverse(), edge1.inverse()].into()],
        surface,
    )
    .unwrap();

    vec![face0, face1, face2, face3].into()
}

#[test]
fn special_cylinder() {
    let shell = special_cylinder_model();
    let mut mesh = shell.triangulation(0.01).to_polygon();
    mesh.put_together_same_attrs(TOLERANCE)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}

#[test]
fn special_cylinder_csolid() {
    let shell = special_cylinder_model().compress();
    let mut mesh = shell.triangulation(0.01).to_polygon();
    mesh.put_together_same_attrs(TOLERANCE)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}

#[test]
fn robust_closed() {
    let cube: Solid = {
        let v = builder::vertex(Point3::origin());
        let e = builder::extrude(&v, Vector3::unit_x());
        let f = builder::extrude(&e, Vector3::unit_y());
        builder::extrude(&f, Vector3::unit_z())
    };

    let o = Point3::new(0.5, 0.5, 0.5);
    cube.edge_iter().for_each(|edge| {
        let curve = edge.curve();

        if let Curve::Line(line) = curve {
            let m = line.subs(0.5);
            let p = m + 0.2 * (o - m);
            let bsp =
                BsplineCurve::new_unchecked(KnotVector::bezier_knot(2), vec![line.0, p, line.1]);
            edge.set_curve(bsp.into());
        }
    });

    assert!(
        cube.triangulation(0.01)
            .face_iter()
            .all(|face| face.surface().is_none())
    );

    let mut mesh = cube.robust_triangulation(0.01).to_polygon();
    mesh.put_together_same_attrs(TOLERANCE2)
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}

/// Verifies that every face mesh from a tessellated solid has consistent
/// winding — all triangles within a face agree on the normal direction.
#[test]
fn per_face_winding_consistency() {
    let shapes: Vec<(&str, Shell)> = vec![
        ("cube", {
            let v = builder::vertex(Point3::origin());
            let e = builder::extrude(&v, Vector3::unit_x());
            let f = builder::extrude(&e, Vector3::unit_y());
            let solid: Solid = builder::extrude(&f, Vector3::unit_z());
            solid.into_boundaries().pop().unwrap()
        }),
        ("cylinder", special_cylinder_model()),
    ];
    for (name, shell) in &shapes {
        let meshed = shell.triangulation(0.01);
        for (fi, face) in meshed.face_iter().enumerate() {
            if let Some(mut poly) = face.surface() {
                if !face.orientation() {
                    poly.invert();
                }
                let positions = poly.positions();
                let normals = poly.normals();
                let (mut pos, mut neg) = (0usize, 0usize);
                for tri in poly.tri_faces() {
                    let e1 = positions[tri[1].pos] - positions[tri[0].pos];
                    let e2 = positions[tri[2].pos] - positions[tri[0].pos];
                    let geo = e1.cross(e2);
                    let vtx: Vector3 = tri.iter().fold(Vector3::zero(), |s, v| {
                        s + v.nor.map(|idx| normals[idx]).unwrap_or(Vector3::zero())
                    });
                    let dot = geo.dot(vtx);
                    if dot > 0.0 {
                        pos += 1;
                    } else if dot < 0.0 {
                        neg += 1;
                    }
                }
                assert!(
                    pos == 0 || neg == 0,
                    "{name} face {fi}: mixed winding ({pos} pos, {neg} neg)"
                );
            }
        }
    }
}
