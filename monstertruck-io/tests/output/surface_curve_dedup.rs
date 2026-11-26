//! Edge associations reference the face's surface entity instead of
//! re-emitting a full copy per edge.
//!
//! The default shell writer associates every edge with the surfaces of its
//! adjacent faces via `SURFACE_CURVE`. Historically each association re-emitted
//! the whole surface entity, so a surface appeared once as the face geometry
//! plus twice per bounding edge (~3x bloat on surface-heavy shells). Each
//! surface must now be emitted exactly once -- as the face geometry -- and the
//! `SURFACE_CURVE` associations must reference that entity.

use monstertruck_io::step::load::Table;
use monstertruck_io::step::save::*;
use monstertruck_modeling::*;

/// Builds a unit cube solid via the extrude chain. Every face is planar, so
/// the emitted surface entity is `PLANE`.
fn unit_cube() -> Solid {
    let v = builder::vertex(Point3::origin());
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    builder::extrude(&f, Vector3::unit_z())
}

/// Counts non-overlapping occurrences of `needle` in `haystack`.
fn count(haystack: &str, needle: &str) -> usize { haystack.matches(needle).count() }

#[test]
fn surface_is_emitted_once_per_face_not_per_edge() {
    let compressed = unit_cube().compress();
    let face_count = compressed.boundaries[0].faces.len();
    assert_eq!(face_count, 6, "a cube has six faces");

    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();

    // The associations still exist (each shared edge lies on two faces).
    assert!(
        step.contains("SURFACE_CURVE"),
        "edges must still associate their surfaces via SURFACE_CURVE:\n{step}"
    );

    // One `PLANE` per face -- no per-edge re-emission. Before the dedup a cube
    // emitted 6 face planes plus two per edge (12 edges) = 30 planes.
    let plane_count = count(&step, "= PLANE(");
    assert_eq!(
        plane_count, face_count,
        "each face surface must be emitted exactly once (got {plane_count} PLANE \
         entities for {face_count} faces):\n{step}"
    );
}

#[test]
fn deduped_surface_curve_still_reloads() {
    let compressed = unit_cube().compress();
    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();

    // Structurally valid and reloads to the same face set.
    assert_eq!(count(&step, "PARAMETRIC_REPRESENTATION_CONTEXT()"), 0);
    let table = Table::from_step(&step).unwrap();
    let step_shell = table
        .shell
        .values()
        .next()
        .expect("reloaded STEP has a shell");
    let reloaded = table.to_compressed_shell(step_shell).unwrap();
    assert_eq!(reloaded.faces.len(), 6, "roundtrip preserves the six faces");
}

#[test]
fn exact_pcurves_reference_their_existing_face_surfaces() {
    let compressed = unit_cube().compress_with_parameter_curves(KERNEL_EPSILON * 100.0);
    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();
    assert_eq!(count(&step, "= PLANE("), 6);
    assert_eq!(count(&step, "= PCURVE("), 24);
    assert_eq!(count(&step, "PARAMETRIC_REPRESENTATION_CONTEXT()"), 1);
    // All parameter representations must refer to the same context.
    // Its identity is shared, not inferred from a coordinate tolerance.
    let context_id = step
        .lines()
        .collect::<Vec<_>>()
        .windows(3)
        .find_map(|lines| {
            (lines[1].trim() == "GEOMETRIC_REPRESENTATION_CONTEXT(2)"
                && lines[2].trim() == "PARAMETRIC_REPRESENTATION_CONTEXT()")
                .then(|| lines[0].split_once(" = ").map(|(id, _)| id))
                .flatten()
        })
        .expect("shared context");
    let context_suffix = format!(", {context_id});");
    assert_eq!(
        step.lines()
            .filter(|line| line.contains("= DEFINITIONAL_REPRESENTATION(")
                && line.ends_with(&context_suffix))
            .count(),
        24
    );
    let table = Table::from_step(&step).expect("deduplicated STEP parses");
    let shell = table.shell.values().next().expect("closed shell");
    let restored = table
        .to_compressed_trimmed_shell(shell)
        .expect("exact trims restore");
    assert_eq!(restored.faces.len(), 6);
    assert_eq!(
        restored
            .faces
            .iter()
            .flat_map(|face| face.boundaries.iter().flatten())
            .filter(|edge| edge.trim_curve.is_some())
            .count(),
        24
    );
}

#[test]
fn shared_context_indices_do_not_overlap_across_models() {
    let first = unit_cube().compress_with_parameter_curves(KERNEL_EPSILON * 100.0);
    let second = unit_cube().compress_with_parameter_curves(KERNEL_EPSILON * 100.0);
    let models: StepModels<_, _, _> = [&first, &second].into_iter().collect();
    let step = CompleteStepDisplay::new(models, Default::default()).to_string();
    assert_eq!(count(&step, "= PCURVE("), 48);
    assert_eq!(count(&step, "PARAMETRIC_REPRESENTATION_CONTEXT()"), 2);
    let entity_ids = step
        .lines()
        .filter_map(|line| line.strip_prefix('#')?.split_once(" = ").map(|(id, _)| id))
        .collect::<Vec<_>>();
    assert_eq!(
        entity_ids
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        entity_ids.len()
    );
    let table = Table::from_step(&step).expect("multiple models parse");
    assert_eq!(table.shell.len(), 2);
    for shell in table.shell.values() {
        let restored = table
            .to_compressed_trimmed_shell(shell)
            .expect("exact trims restore");
        assert_eq!(restored.faces.len(), 6);
        assert_eq!(
            restored
                .faces
                .iter()
                .flat_map(|face| face.boundaries.iter().flatten())
                .filter(|edge| edge.trim_curve.is_some())
                .count(),
            24
        );
    }
}

#[test]
fn pcurve_reference_requires_the_same_emitted_chart() {
    let face = Plane::xy();
    let foreign = Plane::new(
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    );
    let curve = Line(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0));
    let exact = ParameterCurve::new(curve, face);
    let different = ParameterCurve::new(curve, foreign);
    let encoding = StepDisplay::new(face, 1).to_string();
    assert!(exact.matches_surface_encoding(&encoding));
    assert!(!different.matches_surface_encoding(&encoding));
    assert_eq!(
        exact.step_length_with_shared_context(true) + face.step_length() + 1,
        exact.step_length()
    );
    assert_eq!(
        different.step_length_with_shared_context(false) + 1,
        different.step_length()
    );
}
