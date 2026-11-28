use super::*;
use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};

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
    let recovered = recover_proven_test_loft(&step, &kernel)?;
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
