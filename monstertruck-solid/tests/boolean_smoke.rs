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
//! `Solid::try_new`, which only returns `Ok` when every boundary shell is a
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

const TOL: f64 = 0.05;
const STRICT_TOL: f64 = 1.0e-6;
const VOLUME_EPS: f64 = 1.0e-3;

/// Axis-aligned unit cube with its minimum corner at `origin`.
fn unit_cube(origin: Point3) -> Solid {
    let v = builder::vertex(origin);
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    builder::extrude(&f, Vector3::unit_z())
}

/// Square-section vertical column centered on `(cx, cy)`, tall enough to pierce
/// a unit cube standing on `z = 0`.
fn square_column(cx: f64, cy: f64, half: f64) -> Solid {
    let v = builder::vertex(Point3::new(cx - half, cy - half, -0.5));
    let e = builder::extrude(&v, Vector3::unit_x() * (2.0 * half));
    let f = builder::extrude(&e, Vector3::unit_y() * (2.0 * half));
    builder::extrude(&f, Vector3::unit_z() * 2.0)
}

fn unit_cylinder_y(y0: f64) -> Solid {
    let center = Point3::new(0.0, y0, 0.0);
    let vertex = builder::vertex(Point3::new(0.0, y0, 1.0));
    let circle = builder::revolve(
        &vertex,
        center,
        Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    );
    let disk = builder::try_attach_plane(&[circle]).expect("cylinder disk");
    builder::extrude(&disk, Vector3::unit_y())
}

/// Asserts a boolean result is watertight and geometrically correct: it has at
/// least one boundary shell, tessellates to a non-empty mesh, and that mesh has
/// the expected (positively oriented) volume.
fn assert_solid(label: &str, solid: &Solid, expected_volume: f64) -> Result<()> {
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

fn bezier_disk_extrusion(cx: f64, cy: f64, z0: f64, radius: f64, height: f64) -> Result<Solid> {
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

#[test]
fn union_with_thin_overlapping_bezier_protrusion_keeps_host() -> Result<()> {
    let host = unit_cube(Point3::origin());
    let feature = bezier_disk_extrusion(0.5, 0.5, 0.999_99, 0.2, 0.100_01)?;
    let result = monstertruck_solid::or(&host, &feature, STRICT_TOL)?;
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "union must stay consistent"
    );
    let mesh = result.triangulation(0.005).to_polygon();
    let volume = mesh.volume();
    anyhow::ensure!(
        volume > 1.0 && volume < 1.02,
        "union volume {volume:.9} dropped the host or added unrelated material"
    );
    Ok(())
}

#[test]
fn union_of_overlapping_cubes() -> Result<()> {
    // Two unit cubes overlapping in a `0.5` cube: `2 - 0.5^3 = 1.875`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::or(&a, &b, TOL)?;
    assert_solid("union of overlapping cubes", &result, 1.875)
}

#[test]
fn intersection_of_overlapping_cubes() -> Result<()> {
    // The shared `0.5` cube: `0.5^3 = 0.125`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::and(&a, &b, TOL)?;
    assert_solid("intersection of overlapping cubes", &result, 0.125)
}

#[test]
fn difference_of_overlapping_cubes() -> Result<()> {
    // First cube minus the shared region: `1 - 0.125 = 0.875`.
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.5, 0.5, 0.5));
    let result = monstertruck_solid::difference(&a, &b, TOL)?;
    assert_solid("difference of overlapping cubes", &result, 0.875)
}

#[test]
fn cube_minus_column() -> Result<()> {
    // Unit cube with a `0.4 x 0.4` column punched through: `1 - 0.4^2 = 0.84`.
    let cube = unit_cube(Point3::origin());
    let column = square_column(0.5, 0.5, 0.2);
    let result = monstertruck_solid::difference(&cube, &column, STRICT_TOL)?;
    anyhow::ensure!(
        result.boundaries().len() == 1,
        "through-column subtraction must produce one connected boundary shell, got {}",
        result.boundaries().len(),
    );
    anyhow::ensure!(
        result.face_iter().count() == 10,
        "through-column subtraction must produce 10 faces, got {}",
        result.face_iter().count(),
    );
    assert_solid("cube minus square column", &result, 0.84)
}

#[test]
fn profile_generated_difference_matrix() -> Result<()> {
    fn profile_box(min: Point3, max: Point3) -> Result<Solid> {
        let v = builder::vertices([
            Point3::new(min.x, min.y, min.z),
            Point3::new(max.x, min.y, min.z),
            Point3::new(max.x, max.y, min.z),
            Point3::new(min.x, max.y, min.z),
        ]);
        let wire: Wire = vec![
            builder::line(&v[0], &v[1]),
            builder::line(&v[1], &v[2]),
            builder::line(&v[2], &v[3]),
            builder::line(&v[3], &v[0]),
        ]
        .into();
        profile::solid_from_planar_profile::<Curve, Surface>(
            vec![wire],
            Vector3::new(0.0, 0.0, max.z - min.z),
        )
        .map_err(Into::into)
    }

    let primitive_host: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 6.0, 2.0),
    ]));
    let primitive_cutter: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(3.0, 2.0, -1.0),
        Point3::new(7.0, 4.0, 3.0),
    ]));
    let profile_host = profile_box(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 6.0, 2.0))?;
    let profile_cutter = profile_box(Point3::new(3.0, 2.0, -1.0), Point3::new(7.0, 4.0, 3.0))?;
    let profile_cutter_at_origin =
        profile_box(Point3::new(3.0, 2.0, 0.0), Point3::new(7.0, 4.0, 4.0))?;
    let transformed_profile_cutter = builder::transformed(
        &profile_cutter_at_origin,
        Matrix4::from_translation(Vector3::new(0.0, 0.0, -1.0)),
    );

    for (label, solid) in [
        ("primitive_host", &primitive_host),
        ("primitive_cutter", &primitive_cutter),
        ("profile_host", &profile_host),
        ("profile_cutter", &profile_cutter),
        ("transformed_profile_cutter", &transformed_profile_cutter),
    ] {
        eprintln!(
            "input_shell {label} condition={:?} consistent={}",
            solid.boundaries()[0].shell_condition(),
            solid.is_geometric_consistent(),
        );
    }

    for (label, host, cutter) in [
        ("primitive/primitive", &primitive_host, &primitive_cutter),
        ("profile/primitive", &profile_host, &primitive_cutter),
        ("primitive/profile", &primitive_host, &profile_cutter),
        ("profile/profile", &profile_host, &profile_cutter),
        (
            "profile/transformed-profile",
            &profile_host,
            &transformed_profile_cutter,
        ),
    ] {
        match monstertruck_solid::difference(host, cutter, STRICT_TOL) {
            Ok(result) => {
                eprintln!("{label}: ok");
                assert_solid(label, &result, 104.0)?;
            }
            Err(error) => eprintln!("{label}: {error:?}"),
        }
    }

    // Keep the downstream construction as the actual regression assertion.
    let result =
        monstertruck_solid::difference(&profile_host, &transformed_profile_cutter, STRICT_TOL)?;
    assert_solid(
        "profile-generated box minus transformed column",
        &result,
        104.0,
    )
}

#[test]
fn cube_minus_enclosed_cube_preserves_inward_cavity_orientation() -> Result<()> {
    let outer = unit_cube(Point3::origin());
    let inner: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.25, 0.25, 0.25),
        Point3::new(0.75, 0.75, 0.75),
    ]));
    let result = monstertruck_solid::difference(&outer, &inner, TOL)?;
    anyhow::ensure!(
        result.boundaries().len() == 2,
        "enclosed subtraction must produce outer and cavity boundary shells"
    );
    assert_solid("cube minus enclosed cube", &result, 0.875)
}

/// Regression for a radial slot cut from a solid of revolution.
///
/// The cutter reaches the revolution axis and exits through the top end. This
/// is the minimal form of the failure reported by downstream STEP recovery:
/// the classic boolean pipeline used to classify every divided face out of the
/// AND result and return `EmptyOutputShell`.
#[test]
fn revolved_cylinder_minus_axis_touching_radial_slot() -> Result<()> {
    let vertices = builder::vertices([
        Point3::new(0.0, 0.0, -2.0),
        Point3::new(2.0, 0.0, -2.0),
        Point3::new(2.0, 0.0, 2.0),
        Point3::new(0.0, 0.0, 2.0),
    ]);
    // Use the axis-aware wire revolution path. Generic face revolution would
    // sweep the closing on-axis edge into zero-length topology, which is not a
    // valid input for testing the Boolean itself.
    let profile: Wire = vec![
        builder::line(&vertices[0], &vertices[1]),
        builder::line(&vertices[1], &vertices[2]),
        builder::line(&vertices[2], &vertices[3]),
    ]
    .into();
    let shell: Shell = builder::revolve_wire(
        &profile,
        Point3::origin(),
        Vector3::unit_z(),
        builder::SweepAngle::Closed,
        2,
    );
    let host = Solid::try_new(vec![shell])?;
    anyhow::ensure!(
        host.is_geometric_consistent(),
        "axis-aware revolved host is not geometrically consistent"
    );

    let cutter: Solid = primitive::cuboid(BoundingBox::from_iter([
        Point3::new(0.0, -0.25, 1.0),
        Point3::new(2.0001, 0.25, 2.0001),
    ]));

    let result = monstertruck_solid::difference(&host, &cutter, STRICT_TOL)?;
    anyhow::ensure!(
        !result.boundaries().is_empty(),
        "radial-slot difference returned an empty solid"
    );
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "radial-slot difference returned inconsistent topology"
    );
    Ok(())
}

#[test]
fn union_of_end_to_end_cylinders_with_full_circular_interface() -> Result<()> {
    let a = unit_cylinder_y(0.0);
    let b = unit_cylinder_y(1.0);
    let result = monstertruck_solid::or(&a, &b, STRICT_TOL)?;
    anyhow::ensure!(result.boundaries().len() == 1, "union must have one shell");
    anyhow::ensure!(
        result.is_geometric_consistent(),
        "union must be geometrically consistent"
    );
    Ok(())
}

#[test]
fn union_of_coaxial_overlapping_boxes_with_coplanar_sides() -> Result<()> {
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.0, 0.0, 0.5));
    let result = monstertruck_solid::or(&a, &b, 1.0e-6)?;
    assert_solid("union of coaxial overlapping boxes", &result, 1.5)
}

#[test]
fn intersection_of_coaxial_overlapping_boxes_with_coplanar_sides() -> Result<()> {
    let a = unit_cube(Point3::origin());
    let b = unit_cube(Point3::new(0.0, 0.0, 0.5));
    let result = monstertruck_solid::and(&a, &b, 1.0e-6)?;
    assert_solid("intersection of coaxial overlapping boxes", &result, 0.5)
}

#[test]
#[allow(
    clippy::excessive_precision,
    reason = "regression fixture preserves recovered STEP coordinates verbatim"
)]
fn union_with_exact_step_redox_glyph_keeps_host() -> Result<()> {
    let host: Solid = primitive::cuboid(BoundingBox::from_iter([Point3::new(-6.0, -6.0, 0.0), Point3::new(6.0, 6.0, 0.7)]));
    let vertices = builder::vertices([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.067961435629058986, 0.16484251813860773, 0.0),
        Point3::new(0.23424996024462441, 0.23424990508389243, 0.0),
        Point3::new(0.26750765413600686, 0.23135793350139178, 0.0),
        Point3::new(0.331131084124896, 0.29498137038459316, 0.0),
        Point3::new(0.50464954459379419, 0.35137487141653967, 0.0),
        Point3::new(0.70130377772823316, 0.27618353440917343, 0.0),
        Point3::new(0.77938710700311331, 0.16195054655610885, 0.0),
        Point3::new(0.92253980379379108, -0.027473723145671869, 0.0),
        Point3::new(0.86903826745003299, -0.16484250434998149, 0.0),
        Point3::new(0.73745345699281728, -0.22991192702650931, 0.0),
        Point3::new(0.65214024702464313, -0.22991192702650931, 0.0),
        Point3::new(0.61020659701434665, -0.18508634680166747, 0.0),
        Point3::new(0.6550321772378056, -0.14026076657821029, 0.0),
        Point3::new(0.66370807819737454, -0.14170673168444381, 0.0),
        Point3::new(0.7287775008732118, -0.14170673168444381, 0.0),
        Point3::new(0.83288858818614253, -0.028919743411946364, 0.0),
        Point3::new(0.72010159991364375, 0.083867244860552859, 0.0),
        Point3::new(0.7085337135808718, 0.083867244860552859, 0.0),
        Point3::new(0.50609550970072226, 0.26461564118209147, 0.0),
        Point3::new(0.31233326193948585, 0.12580095003019576, 0.0),
        Point3::new(0.23569595293087886, 0.14749075758950081, 0.0),
        Point3::new(0.091097235873924376, 0.0028920405325472132, 0.0),
        Point3::new(0.1778564661083748, -0.13013879019162999, 0.0),
        Point3::new(0.31377925462643397, -0.034703631418292247, 0.0),
        Point3::new(0.42222829241949533, -0.083867189701898326, 0.0),
        Point3::new(0.615990567760059, 0.017351912237826284, 0.0),
        Point3::new(0.66081614798420985, -0.059285424350106197, 0.0),
        Point3::new(0.45837794410406074, -0.16628846945622211, 0.0),
        Point3::new(0.45837794410406074, -0.18074834116288763, 0.0),
        Point3::new(0.31377922704641437, -0.32534705821983589, 0.0),
        Point3::new(0.17641043205139795, -0.22701994165262107, 0.0),
        Point3::new(0.060731458406389827, -0.15761256849734195, 0.0),
    ]);
    let wire: Wire = vec![
        builder::bezier(&vertices[0], &vertices[1], vec![Point3::new(0.0, 0.062177437301247807, 0.0), Point3::new(0.02458178258896071, 0.12146291681139321, 0.0)]),
        builder::bezier(&vertices[1], &vertices[2], vec![Point3::new(0.1127870434325362, 0.20966812594347317, 0.0), Point3::new(0.17207250915197969, 0.23424990508389243, 0.0)]),
        builder::bezier(&vertices[2], &vertices[3], vec![Point3::new(0.24581784657739414, 0.23424990508389243, 0.0), Point3::new(0.25593976780323668, 0.23280392618764623, 0.0)]),
        builder::bezier(&vertices[3], &vertices[4], vec![Point3::new(0.28630547632141701, 0.25593971264319659, 0.0), Point3::new(0.30654929119308161, 0.27762952020111697, 0.0)]),
        builder::bezier(&vertices[4], &vertices[5], vec![Point3::new(0.38174063509545153, 0.33113105137345045, 0.0), Point3::new(0.44102611460490326, 0.35137487141653967, 0.0)]),
        builder::bezier(&vertices[5], &vertices[6], vec![Point3::new(0.57694890312227143, 0.35137487141653967, 0.0), Point3::new(0.64635627627823844, 0.32390111552080825, 0.0)]),
        builder::bezier(&vertices[6], &vertices[7], vec![Point3::new(0.7360074918858901, 0.24437181252041817, 0.0), Point3::new(0.76348127019090706, 0.20533016167194695, 0.0)]),
        builder::bezier(&vertices[7], &vertices[8], vec![Point3::new(0.86180833159739301, 0.13592278851667228, 0.0), Point3::new(0.92109383868686567, 0.059285451928738908, 0.0)]),
        builder::bezier(&vertices[8], &vertices[9], vec![Point3::new(0.92253980379379108, -0.078083274115535062, 0.0), Point3::new(0.90374198160838404, -0.12724685997916563, 0.0)]),
        builder::bezier(&vertices[9], &vertices[10], vec![Point3::new(0.83433460845310536, -0.20243814872080623, 0.0), Point3::new(0.78806300796198947, -0.22557397654638134, 0.0)]),
        builder::line(&vertices[10], &vertices[11]),
        builder::bezier(&vertices[11], &vertices[12], vec![Point3::new(0.62900441919906358, -0.22846596192027668, 0.0), Point3::new(0.61020659701434665, -0.20822211946720426, 0.0)]),
        builder::bezier(&vertices[12], &vertices[13], vec![Point3::new(0.61020659701434665, -0.16050460902989006, 0.0), Point3::new(0.6304503843059921, -0.14026076657821029, 0.0)]),
        builder::bezier(&vertices[13], &vertices[14], vec![Point3::new(0.6579241074509703, -0.14026076657821029, 0.0), Point3::new(0.66081609282417242, -0.14170673168444381, 0.0)]),
        builder::line(&vertices[14], &vertices[15]),
        builder::bezier(&vertices[15], &vertices[16], vec![Point3::new(0.78806300796198947, -0.13736878120431495, 0.0), Point3::new(0.83288858818614253, -0.088205222922089099, 0.0)]),
        builder::bezier(&vertices[16], &vertices[17], vec![Point3::new(0.83288858818614253, 0.033257693890688778, 0.0), Point3::new(0.78227903721627801, 0.083867244860552859, 0.0)]),
        builder::bezier(&vertices[17], &vertices[18], vec![Point3::new(0.71576364943351312, 0.083867244860552859, 0.0), Point3::new(0.71287166406030922, 0.083867244860552859, 0.0)]),
        builder::bezier(&vertices[18], &vertices[19], vec![Point3::new(0.69696582724810296, 0.18508634680028102, 0.0), Point3::new(0.61020659701434665, 0.26461564118209147, 0.0)]),
        builder::bezier(&vertices[19], &vertices[20], vec![Point3::new(0.41499832898614297, 0.26461564118209147, 0.0), Point3::new(0.33980698508446405, 0.20677615435958696, 0.0)]),
        builder::bezier(&vertices[20], &vertices[21], vec![Point3::new(0.2906434543815668, 0.13881482905060594, 0.0), Point3::new(0.263169703656569, 0.14749075758950081, 0.0)]),
        builder::bezier(&vertices[21], &vertices[22], vec![Point3::new(0.15616665854976075, 0.14749075758950081, 0.0), Point3::new(0.091097235873924376, 0.08242133491435677, 0.0)]),
        builder::bezier(&vertices[22], &vertices[23], vec![Point3::new(0.091097235873924376, -0.056393438977597299, 0.0), Point3::new(0.12724691513851072, -0.10700298994746138, 0.0)]),
        builder::bezier(&vertices[23], &vertices[24], vec![Point3::new(0.19810028098073218, -0.07519126116162056, 0.0), Point3::new(0.25160178974377923, -0.034703631418292247, 0.0)]),
        builder::bezier(&vertices[24], &vertices[25], vec![Point3::new(0.35715886974365851, -0.034703631418292247, 0.0), Point3::new(0.39620053438075198, -0.053501453603701954, 0.0)]),
        builder::line(&vertices[25], &vertices[26]),
        builder::line(&vertices[26], &vertices[27]),
        builder::line(&vertices[27], &vertices[28]),
        builder::bezier(&vertices[28], &vertices[29], vec![Point3::new(0.45837794410406074, -0.17207244020262191, 0.0), Point3::new(0.45837794410406074, -0.17641039068275077, 0.0)]),
        builder::bezier(&vertices[29], &vertices[30], vec![Point3::new(0.45837794410406074, -0.26027763554469363, 0.0), Point3::new(0.3933085214275307, -0.32534705821983589, 0.0)]),
        builder::bezier(&vertices[30], &vertices[31], vec![Point3::new(0.25015576947681151, -0.32534705821983589, 0.0), Point3::new(0.19520826802681723, -0.28485942847650936, 0.0)]),
        builder::bezier(&vertices[31], &vertices[32], vec![Point3::new(0.13158483803792853, -0.21545205531985179, 0.0), Point3::new(0.092543173400141399, -0.19231628265569256, 0.0)]),
        builder::bezier(&vertices[32], &vertices[0], vec![Point3::new(0.021689807557919494, -0.11423300854015572, 0.0), Point3::new(0.0, -0.057839486823890951, 0.0)]),
    ].into();
    let feature = profile::solid_from_planar_profile::<Curve, Surface>(vec![wire], Vector3::new(0.0, 0.0, 0.010010000000000008))?;
    let feature = builder::transformed(&feature, Matrix4::from_translation(Vector3::new(2.8117697409035993, 4.5510814676812847, 0.69999)));
    let result = monstertruck_solid::or(&host, &feature, STRICT_TOL)?;
    anyhow::ensure!(result.is_geometric_consistent(), "exact glyph union must stay consistent");
    let volume = result.triangulation(0.005).to_polygon().volume();
    eprintln!("exact glyph union faces={} volume={volume:.9}", result.face_iter().count());
    anyhow::ensure!(volume > 100.8 && volume < 101.0, "exact glyph union volume {volume:.9} dropped the host or added unrelated material");
    Ok(())
}
