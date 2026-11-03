//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;
use crate::step::load::step_geometry::SurfaceCurveRepresentation as StepSurfaceCurveRepresentation;
use crate::step::load::step_geometry::{
    Curve2D, ElementarySurface, SurfaceCurve3D, SurfaceCurveAssociatedGeometry,
};
use monstertruck_geometry::prelude::{
    BoundedCurve, Line, ParametricCurve, ParametricSurface, Processor, RevolutionSurface,
    Tolerance, Vector3,
};
use std::f64::consts::TAU;

/// Loads `occt-cylinder.step` and builds its trimmed shell. The
/// cylinder fixture carries `PCURVE` entities on both planar and
/// cylindrical surfaces, so this exercises the `ToSameGeometry<Curve2D>`
/// load path end-to-end -- the unit tests in `step_geometry/geom_impls`
/// only verify the impls in isolation. At least one of the resulting
/// trim curves must be present (`Some(_)`) and at least one of them
/// must contain a 2D curve variant that comes through the conversion.
#[test]
fn pcurve_load_path_populates_trim_curves() -> anyhow::Result<()> {
    let step_string = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../resources/step/occt-cylinder.step",
    ));
    let table = crate::step::load::Table::from_step(step_string)?;
    let step_shell = table
        .shell
        .values()
        .next()
        .ok_or_else(|| anyhow::anyhow!("the cylinder fixture must contain a STEP shell."))?;
    let trimmed = table.to_compressed_trimmed_shell(step_shell)?;
    let total_edge_uses: usize = trimmed
        .faces
        .iter()
        .flat_map(|face| &face.boundaries)
        .map(|wire| wire.len())
        .sum();
    let trim_curves_present: usize = trimmed
        .faces
        .iter()
        .flat_map(|face| &face.boundaries)
        .flat_map(|wire| wire.iter())
        .filter(|edge_use| edge_use.trim_curve.is_some())
        .count();
    assert!(
        total_edge_uses > 0,
        "the cylinder fixture should have at least one edge-use after trimmed loading.",
    );
    assert!(
        trim_curves_present > 0,
        "at least one edge-use should carry a trim curve. \
         total edge-uses: {total_edge_uses}.",
    );
    Ok(())
}

#[test]
fn trimmed_solid_report_is_aligned_with_boundaries() -> anyhow::Result<()> {
    let step_string = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../resources/step/occt-cylinder.step",
    ));
    let table = crate::step::load::Table::from_step(step_string)?;
    let step_solid = table
        .manifold_solid_brep
        .values()
        .next()
        .ok_or_else(|| anyhow::anyhow!("the cylinder fixture must contain a STEP solid."))?;

    let (trimmed, reports) = table.to_compressed_trimmed_solid_reported(step_solid)?;
    assert_eq!(trimmed.boundaries.len(), reports.len());
    assert_eq!(trimmed.boundaries.len(), 1);
    assert_eq!(
        table.to_compressed_trimmed_solid(step_solid)?,
        trimmed,
        "the compatibility wrapper must discard only the report",
    );
    Ok(())
}

#[test]
fn trimmed_solid_import_assigns_stable_identity_and_provenance() -> anyhow::Result<()> {
    let step_string = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../resources/step/occt-cylinder.step",
    ));
    let table = crate::step::load::Table::from_step(step_string)?;
    let step_solid = table
        .manifold_solid_brep
        .values()
        .next()
        .ok_or_else(|| anyhow::anyhow!("the cylinder fixture must contain a STEP solid."))?;

    let import = table.to_compressed_trimmed_solid_import(step_solid)?;
    assert_eq!(import.solid.boundaries.len(), import.shell_reports.len());

    let mut vertex_count = 0usize;
    let mut edge_count = 0usize;
    let mut face_count = 0usize;
    let mut max_stable_id = 0u64;
    for shell in &import.solid.boundaries {
        let vertex_ids = shell
            .vertex_stable_ids
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("imported vertices must carry stable IDs"))?;
        let edge_ids = shell
            .edge_stable_ids
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("imported edges must carry stable IDs"))?;
        let face_ids = shell
            .face_stable_ids
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("imported faces must carry stable IDs"))?;
        assert_eq!(vertex_ids.len(), shell.vertices.len());
        assert_eq!(edge_ids.len(), shell.edges.len());
        assert_eq!(face_ids.len(), shell.faces.len());
        for id in vertex_ids.iter().chain(edge_ids).chain(face_ids) {
            assert!(id.is_assigned());
            max_stable_id = max_stable_id.max(id.raw());
        }
        vertex_count += shell.vertices.len();
        edge_count += shell.edges.len();
        face_count += shell.faces.len();
    }

    assert_eq!(import.provenance.vertices.len(), vertex_count);
    assert_eq!(import.provenance.edges.len(), edge_count);
    assert_eq!(import.provenance.faces.len(), face_count);
    assert!(
        import
            .provenance
            .vertices
            .iter()
            .all(|entry| entry.stable_id.is_assigned())
    );
    assert!(
        import
            .provenance
            .edges
            .iter()
            .all(|entry| entry.stable_id.is_assigned())
    );
    assert!(
        import
            .provenance
            .faces
            .iter()
            .all(|entry| entry.stable_id.is_assigned())
    );
    assert!(
        import
            .solid
            .id_allocator
            .as_ref()
            .is_some_and(|allocator| allocator.peek() > max_stable_id),
        "the retained allocator must continue above every imported stable ID",
    );
    Ok(())
}

#[test]
fn plain_and_trimmed_solid_imports_assign_the_same_source_identity() -> anyhow::Result<()> {
    let step_string = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../resources/step/occt-cylinder.step",
    ));
    let table = crate::step::load::Table::from_step(step_string)?;
    let step_solid = table
        .manifold_solid_brep
        .values()
        .next()
        .ok_or_else(|| anyhow::anyhow!("the cylinder fixture must contain a STEP solid."))?;

    let plain = table.to_compressed_solid_import(step_solid)?;
    let trimmed = table.to_compressed_trimmed_solid_import(step_solid)?;
    assert_eq!(plain.provenance, trimmed.provenance);
    assert_eq!(plain.shell_reports, trimmed.shell_reports);
    Ok(())
}

#[test]
fn stable_id_assigner_reuses_source_identity_only_for_named_sources() -> anyhow::Result<()> {
    let mut assigner = StableIdAssigner::default();
    let named = assigner.id_for_source(Some(42))?;
    assert_eq!(assigner.id_for_source(Some(42))?, named);

    let first_owned = assigner.id_for_source(None)?;
    let second_owned = assigner.id_for_source(None)?;
    assert_ne!(named, first_owned);
    assert_ne!(first_owned, second_owned);
    Ok(())
}

fn cylinder_surface() -> Surface {
    let axis = Vector3::unit_z();
    let center = Point3::origin();
    let point = Point3::new(1.0, 0.0, 0.0);
    let line = Line(point, point + axis);
    Surface::Elementary(ElementarySurface::CylindricalSurface(Processor::new(
        RevolutionSurface::by_revolution(line, center, axis),
    )))
}

fn line_pcurve(surface: &Surface, u: f64) -> step_geometry::StepParameterCurve {
    step_geometry::StepParameterCurve::new(
        Box::new(Curve2D::Line(Line(
            Point2::new(u, 0.0),
            Point2::new(u, 1.0),
        ))),
        Box::new(surface.clone()),
    )
}

fn seam_curve(surface: &Surface) -> Curve3D {
    let leader = Curve3D::Line(Line(surface.evaluate(0.0, 0.0), surface.evaluate(0.0, 1.0)));
    Curve3D::SurfaceCurve(SurfaceCurve3D::new(
        StepSurfaceCurveKind::Seam,
        Box::new(leader),
        vec![
            SurfaceCurveAssociatedGeometry::ParameterCurve(line_pcurve(surface, 0.0)),
            SurfaceCurveAssociatedGeometry::ParameterCurve(line_pcurve(surface, TAU)),
        ],
        StepSurfaceCurveRepresentation::ParameterCurve1,
    ))
}

#[test]
fn seam_curve_opposite_orientations_use_opposite_parameter_curves() {
    let surface = cylinder_surface();
    let curve = seam_curve(&surface);

    let forward = Table::exact_trim_curve_on(&curve, &surface, true)
        .expect("forward trim curve should exist");
    let backward = Table::exact_trim_curve_on(&curve, &surface, false)
        .expect("backward trim curve should exist");

    let forward_start = forward.curve().evaluate(forward.curve().range_tuple().0);
    let backward_start = backward.curve().evaluate(backward.curve().range_tuple().0);

    assert!(forward_start.x.near(&0.0));
    assert!(backward_start.x.near(&TAU));
}
