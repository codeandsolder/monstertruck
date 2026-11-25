//! A boundary-constrained recovery patch must survive normal tessellation,
//! including its collapsed parameter edge and shared radial topology.

use monstertruck_geometry::prelude::{BsplineCurve, KnotVector, NurbsCurve};
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::{
    Curve, Edge, Point3, Shell, Vector4, builder, mesh_recovery::approximate_boundary_facet,
};

#[test]
fn reconstructed_curved_facet_tessellates_without_dropped_faces() -> anyhow::Result<()> {
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    let [a, b, c] = points.map(builder::vertex);
    let curve = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
        KnotVector::bezier_knot(2),
        vec![
            Vector4::new(0.0, 0.0, 0.0, 1.0),
            Vector4::new(0.7, 0.0, 0.14, 0.7),
            Vector4::new(2.0, 0.0, 0.0, 1.0),
        ],
    )?));
    let edges = [
        Edge::new(&a, &b, curve)?,
        builder::line(&b, &c),
        builder::line(&c, &a),
    ];
    let approximation = approximate_boundary_facet(&points, &edges, 0.2)?;
    let shell = approximation.faces.into_iter().collect::<Shell>();
    let mesh = shell.triangulation(0.01);
    assert_eq!(mesh.face_iter().count(), 3);
    for face in mesh.face_iter() {
        let polygon = face.surface().expect("no dropped face");
        assert!(!polygon.tri_faces().is_empty());
        assert!(polygon.positions().iter().all(|p| p.z.abs() <= 0.2));
    }
    Ok(())
}

#[test]
fn closed_hybrid_shell_preserves_analytic_carrier_and_shared_edges() -> anyhow::Result<()> {
    use monstertruck_modeling::{
        Face, Matrix4, Processor, Sphere, Surface, Vector3, Wire,
        mesh_recovery::{
            BoundaryPatchOptions, PrescribedSourceEdge, SurfaceCarrier, SurfaceRegionBudget,
            approximate_boundary_patch,
        },
    };

    // Three orthogonal sphere points in a compact longitude/latitude chart.
    // Rotation avoids both poles and makes this an ordinary analytic patch,
    // not a special spherical-cap chart.
    use monstertruck_geometry::prelude::{InnerSpace, SquareMatrix, Transform};
    let axis = Vector3::new(1.0, -1.0, 0.0).normalize();
    let tangent = Vector3::new(1.0, 1.0, -2.0).normalize();
    let normal = tangent.cross(axis);
    let matrix = Matrix4::from_cols(
        normal.extend(0.0),
        tangent.extend(0.0),
        axis.extend(0.0),
        Vector4::new(0.0, 0.0, 0.0, 1.0),
    );
    let inverse = matrix
        .invert()
        .ok_or_else(|| anyhow::anyhow!("chart frame"))?;
    let source = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, 0.0),
    ]
    .map(|p| inverse.transform_point(p));
    let carrier = SurfaceCarrier::CircularRevolution {
        axial_center: Point3::new(0.0, 0.0, 0.0),
        axis: Vector3::unit_z(),
        radial_center: 0.0,
        profile_radius: 1.0,
    };
    let budget = SurfaceRegionBudget {
        vertex: 1.0e-10,
        chord: 0.35,
    };
    let vertices = source.map(builder::vertex);
    let mut prescribed = Vec::new();
    for [a, b] in [[0, 1], [1, 2], [2, 0]] {
        let lifted = carrier
            .lift_source_chord([source[a], source[b]], budget)
            .ok_or_else(|| anyhow::anyhow!("sphere chord"))?;
        let edge = Edge::new(&vertices[a], &vertices[b], Curve::NurbsCurve(lifted.curve))?;
        prescribed.push(PrescribedSourceEdge {
            vertices: [a, b],
            edge,
        });
    }
    let analytic = Face::new(
        vec![prescribed.iter().map(|p| p.edge.clone()).collect::<Wire>()],
        Surface::Sphere(Processor::new(Sphere::new(Point3::new(0.0, 0.0, 0.0), 1.0))),
    )?;
    let fallback_polygons: [&[usize]; 3] = [&[1, 0, 3], &[2, 1, 3], &[0, 2, 3]];
    let approximation = approximate_boundary_patch(
        &source,
        &fallback_polygons,
        &prescribed,
        BoundaryPatchOptions::new(budget),
    )?;
    let mut shell = approximation.shell;
    shell.push(analytic);
    shell.check_solid_boundary()?;
    assert!(shell.is_geometric_consistent());
    let mesh = shell.triangulation(0.005);
    assert_eq!(mesh.face_iter().count(), 10);
    for face in mesh.face_iter() {
        let polygon = face.surface().expect("hybrid face tessellates");
        assert!(!polygon.tri_faces().is_empty());
    }
    let mut polygon = mesh.to_polygon();
    polygon
        .put_together_same_attrs(monstertruck_modeling::KERNEL_EPSILON)
        .remove_degenerate_faces();
    assert_eq!(
        polygon.shell_condition(),
        monstertruck_topology::shell::ShellCondition::Closed
    );
    // The three radial sides are planar and the fourth face is the
    // unit sphere octant. Its positive volume is exactly pi/6.
    let expected = std::f64::consts::PI / 6.0;
    assert!((polygon.volume() - expected).abs() < expected * 0.03);
    Ok(())
}
