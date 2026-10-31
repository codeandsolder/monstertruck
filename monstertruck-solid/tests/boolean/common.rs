//! Downstream-application smoke scenario: build primitive solids, combine them
//! with the boolean operators, and confirm the results tessellate into
//! non-degenerate watertight meshes.
//!
//! This mirrors how an interactive CAD app (e.g. an OpenCADStudio-style scene
//! graph) drives the kernel: a handful of primitives fed through
//! `union`/`intersection`/`difference`, then meshed for display. It is a
//! coarse "does the pipeline stay closed end to end" check, not a numerical
//! accuracy test -- the `transversal::integrate` unit tests own the exact
//! values.
//!
//! Note: `monstertruck_solid::{and, or}` build their output through
//! `Solid::new`, which only returns `Ok` when every boundary shell is a
//! closed manifold. So a returned `Solid` is watertight by construction; these
//! tests additionally confirm it tessellates to a non-empty mesh whose volume
//! matches the analytic result. Because every primitive here is flat-faced, the
//! tessellated volume is exact (independent of mesh density), so an exact-value
//! assertion is both legitimate and a far stronger regression guard than a
//! "non-empty output" smoke check -- it would catch a boolean that silently
//! kept the wrong region.

use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::*;

pub(super) const TOL: f64 = 0.05;
pub(super) const STRICT_TOL: f64 = 1.0e-6;
pub(super) const VOLUME_EPS: f64 = 1.0e-3;

/// Axis-aligned unit cube with its minimum corner at `origin`.
pub(super) fn unit_cube(origin: Point3) -> Solid {
    let v = builder::vertex(origin);
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    builder::extrude(&f, Vector3::unit_z())
}

/// Square-section vertical column centered on `(cx, cy)`, tall enough to pierce
/// a unit cube standing on `z = 0`.
pub(super) fn square_column(cx: f64, cy: f64, half: f64) -> Solid {
    let v = builder::vertex(Point3::new(cx - half, cy - half, -0.5));
    let e = builder::extrude(&v, Vector3::unit_x() * (2.0 * half));
    let f = builder::extrude(&e, Vector3::unit_y() * (2.0 * half));
    builder::extrude(&f, Vector3::unit_z() * 2.0)
}

pub(super) fn unit_cylinder_y(y0: f64) -> Result<Solid> {
    let center = Point3::new(0.0, y0, 0.0);
    let vertex = builder::vertex(Point3::new(0.0, y0, 1.0));
    let circle = builder::revolve(
        &vertex,
        center,
        Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    )?;
    let disk = builder::try_attach_plane([circle])?;
    Ok(builder::extrude(&disk, Vector3::unit_y()))
}

/// Asserts a boolean result is watertight and geometrically correct: it has at
/// least one boundary shell, tessellates to a non-empty mesh, and that mesh has
/// the expected (positively oriented) volume.
pub(super) fn assert_solid(label: &str, solid: &Solid, expected_volume: f64) -> Result<()> {
    anyhow::ensure!(
        !solid.boundaries().is_empty(),
        "{label}: solid has no boundary shells"
    );
    let mesh = solid.triangulation(0.01).to_polygon();
    let triangles = mesh.faces().triangle_iter().count();
    anyhow::ensure!(triangles > 0, "{label}: tessellation produced no triangles");
    let volume = mesh.volume();
    anyhow::ensure!(
        (volume - expected_volume).abs() < VOLUME_EPS,
        "{label}: tessellated volume {volume:.6}, expected {expected_volume:.6}"
    );
    Ok(())
}

/// Independent point-in-polyhedron oracle for Boolean tests.
///
/// Production Boolean fallbacks use signed ray crossings through
/// `PolygonMesh::inside`. Reusing that here would let the test and the code
/// under test share the same classifier defect. Instead sum the oriented
/// solid angle of every triangle: a closed oriented mesh winds by ±4π around
/// an interior point and approximately zero around an exterior point.
fn solid_angle_contains(mesh: &PolygonMesh, point: Point3) -> bool {
    let mut winding = 0.0;
    for triangle in mesh.faces().triangle_iter() {
        let a = mesh.positions()[triangle[0].pos] - point;
        let b = mesh.positions()[triangle[1].pos] - point;
        let c = mesh.positions()[triangle[2].pos] - point;
        let la = a.magnitude();
        let lb = b.magnitude();
        let lc = c.magnitude();
        let numerator = a.dot(b.cross(c));
        let denominator = la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
        winding += 2.0 * numerator.atan2(denominator);
    }
    winding.abs() > std::f64::consts::TAU
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct UnionOccupancyStats {
    pub lhs_only: usize,
    pub rhs_only: usize,
    pub overlap: usize,
}

/// Checks a union as a set operation instead of pinning one particular B-rep
/// decomposition or a reference kernel's answer.
///
/// The two operands are meshed *before* the Boolean, and all three meshes are
/// classified with the independent solid-angle oracle above rather than the
/// production ray-crossing classifier. The
/// result must agree with `lhs || rhs` over deterministic samples spanning
/// each operand and the combined extent. Volume is only used for the universal
/// set bounds `max(Va, Vb) <= Vunion <= Va + Vb`, never as a captured golden.
pub(super) fn assert_union_matches_operand_occupancy(
    label: &str,
    lhs: &Solid,
    rhs: &Solid,
    union: &Solid,
) -> Result<UnionOccupancyStats> {
    let lhs_mesh = lhs.triangulation(0.005).to_polygon();
    let rhs_mesh = rhs.triangulation(0.005).to_polygon();
    let union_mesh = union.triangulation(0.005).to_polygon();

    let lhs_volume = lhs_mesh.volume().abs();
    let rhs_volume = rhs_mesh.volume().abs();
    let union_volume = union_mesh.volume().abs();
    let volume_slack = 1.0e-3 * (1.0 + lhs_volume + rhs_volume);
    anyhow::ensure!(
        union_volume + volume_slack >= lhs_volume.max(rhs_volume),
        "{label}: union volume {union_volume:.9} is smaller than an operand ({lhs_volume:.9}, {rhs_volume:.9})"
    );
    anyhow::ensure!(
        union_volume <= lhs_volume + rhs_volume + volume_slack,
        "{label}: union volume {union_volume:.9} exceeds the sum of operand volumes ({lhs_volume:.9} + {rhs_volume:.9})"
    );

    let lhs_box = lhs_mesh.bounding_box();
    let rhs_box = rhs_mesh.bounding_box();
    let mut combined_box = lhs_box;
    combined_box += rhs_box;
    let boxes = [lhs_box, rhs_box, combined_box];
    // Mid-cell-ish, deliberately asymmetric fractions avoid the many exact
    // planes/axes in CAD fixtures without introducing randomness.
    let fractions = [0.113, 0.271, 0.419, 0.587, 0.743, 0.901];

    let mut stats = UnionOccupancyStats::default();
    let mut checked = 0usize;
    for bbox in boxes {
        let min = bbox.min();
        let diagonal = bbox.diagonal();
        for fx in fractions {
            for fy in fractions {
                for fz in fractions {
                    let point = Point3::new(
                        min.x + diagonal.x * fx,
                        min.y + diagonal.y * fy,
                        min.z + diagonal.z * fz,
                    );
                    let in_lhs = solid_angle_contains(&lhs_mesh, point);
                    let in_rhs = solid_angle_contains(&rhs_mesh, point);
                    let expected = in_lhs || in_rhs;
                    let actual = solid_angle_contains(&union_mesh, point);
                    anyhow::ensure!(
                        actual == expected,
                        "{label}: union occupancy disagrees at {point:?}: lhs={in_lhs}, rhs={in_rhs}, union={actual}"
                    );
                    match (in_lhs, in_rhs) {
                        (true, false) => stats.lhs_only += 1,
                        (false, true) => stats.rhs_only += 1,
                        (true, true) => stats.overlap += 1,
                        (false, false) => {}
                    }
                    checked += 1;
                }
            }
        }
    }
    anyhow::ensure!(checked > 0, "{label}: occupancy oracle sampled no points");
    Ok(stats)
}

pub(super) fn bezier_disk_extrusion(
    cx: f64,
    cy: f64,
    z0: f64,
    radius: f64,
    height: f64,
) -> Result<Solid> {
    let k = 0.552_284_749_830_793_6 * radius;
    let p0 = builder::vertex(Point3::new(cx + radius, cy, z0));
    let p1 = builder::vertex(Point3::new(cx, cy + radius, z0));
    let p2 = builder::vertex(Point3::new(cx - radius, cy, z0));
    let p3 = builder::vertex(Point3::new(cx, cy - radius, z0));
    let wire: Wire = vec![
        builder::bezier(
            &p0,
            &p1,
            vec![
                Point3::new(cx + radius, cy + k, z0),
                Point3::new(cx + k, cy + radius, z0),
            ],
        ),
        builder::bezier(
            &p1,
            &p2,
            vec![
                Point3::new(cx - k, cy + radius, z0),
                Point3::new(cx - radius, cy + k, z0),
            ],
        ),
        builder::bezier(
            &p2,
            &p3,
            vec![
                Point3::new(cx - radius, cy - k, z0),
                Point3::new(cx - k, cy - radius, z0),
            ],
        ),
        builder::bezier(
            &p3,
            &p0,
            vec![
                Point3::new(cx + k, cy - radius, z0),
                Point3::new(cx + radius, cy - k, z0),
            ],
        ),
    ]
    .into();
    profile::solid_from_planar_profile::<Curve, Surface>(vec![wire], Vector3::new(0.0, 0.0, height))
        .map_err(Into::into)
}
