use super::*;
use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
use crate::step_graph::build_index;
use monstertruck_io::step::load::ParsedStepSource;

fn assert_polygon_loft(points: &[[f64; 2]]) -> anyhow::Result<()> {
    let kernel = MonstertruckKernel::project_standard();
    let mut model = CadModel::new();
    let sections = [1.0, 0.9, 0.8]
        .into_iter()
        .zip([0.0, 1.0, 3.0])
        .map(|(scale, level)| {
            Ok(LoftSection {
                profile: Profile2d::polygon(
                    points.iter().map(|[x, y]| [x * scale, y * scale]).collect(),
                )?,
                transform: RigidTransform::translation_mm([17.0, -23.0, level]),
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let root = model.add_node(CadNode::Loft { sections });
    model.add_root(root)?;
    let authored = kernel.evaluate(&model, root)?;
    let step = kernel.to_step(&authored)?;
    let parsed = ParsedStepSource::parse(step.as_bytes())?.exchange;
    let entities = &parsed.data[0].entities;
    let index = build_index(entities);
    let ids = crate::brep::manifold_solid_ids(entities);
    let [id] = ids.as_slice() else {
        anyhow::bail!("expected a single source solid");
    };
    let reduced = ReducedSolid::from_step(*id, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("source B-rep missing"))?;
    let recovered = detect_sectioned_planar_loft(&reduced)
        .ok_or_else(|| anyhow::anyhow!("source ring correspondence not recognized"))?;
    let CadNode::Loft { sections } = recovered.model.node(recovered.root)? else {
        anyhow::bail!("expected an editable loft");
    };
    assert_eq!(sections.len(), 3);
    assert!(sections.iter().all(|section| {
        section
            .profile
            .single_polygon_points()
            .is_some_and(|p| p.len() == points.len())
    }));
    let proof = crate::cad_recovery::SourceRecovery::new(*id, entities, &index, &kernel)
        .ok_or_else(|| anyhow::anyhow!("source boundary certificate unavailable"))?;
    assert!(proof.proves(&recovered.model, recovered.root, None)?);
    Ok(())
}

#[test]
fn triangular_and_nonorthogonal_pentagonal_rings() -> anyhow::Result<()> {
    assert_polygon_loft(&[[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]])?;
    assert_polygon_loft(&[[0.0, 0.0], [3.0, 0.0], [4.0, 2.0], [2.0, 4.0], [-1.0, 2.0]])
}

#[test]
fn concave_rings_keep_source_correspondence() -> anyhow::Result<()> {
    assert_polygon_loft(&[
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 1.0],
        [1.0, 1.0],
        [1.0, 3.0],
        [0.0, 3.0],
    ])
}
