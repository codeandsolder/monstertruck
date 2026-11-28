use super::*;

#[test]
fn cycle_order_is_deterministic() -> anyhow::Result<()> {
    let cylinder = Carrier::Cylinder(CylinderRun {
        origin_mm: [0.0; 3],
        direction: [1.0, 0.0, 0.0],
        radius_mm: 0.1,
        face_ids: vec![],
    });
    let circle = CircleSupport {
        center_mm: [0.0; 3],
        normal: [1.0, 0.0, 0.0],
        x_direction: [0.0, 1.0, 0.0],
        radius_mm: 0.1,
    };
    let graph = Graph {
        pieces: vec![cylinder; 4],
        interfaces: (0..4)
            .map(|i| Interface {
                groups: [i, (i + 1) % 4],
                circle,
            })
            .collect(),
        adjacency: (0..4).map(|i| vec![(i + 3) % 4, i]).collect(),
        radius: 0.1,
    };
    let path = graph
        .walk(3, 0, true)
        .ok_or_else(|| anyhow::anyhow!("cycle missing"))?;
    assert_eq!(
        path.iter().map(|p| p.0).collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );
    assert_eq!(
        path.iter().map(|p| p.2).collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );
    Ok(())
}

#[cfg(feature = "cad-kernel-monstertruck")]
#[test]
fn triangular_spine_with_independent_radii_roundtrips() -> anyhow::Result<()> {
    use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
    use crate::cad_recovery::{SourceRecovery, recover_round_sweep_fragment};
    use crate::math3::add;
    use crate::step_graph::build_index;
    use monstertruck_io::step::load::ParsedStepSource;
    let vertices = [[0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [3.0, 6.0, 0.0]];
    let axis = [0.0, 0.0, 1.0];
    let mut corners = Vec::new();
    for (i, radius) in [0.4, 0.7, 0.9].into_iter().enumerate() {
        let v = vertices[i];
        let incoming = normalize(sub(v, vertices[(i + 2) % 3]), GEOM_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("incoming tangent missing"))?;
        let outgoing = normalize(sub(vertices[(i + 1) % 3], v), GEOM_TOL_MM)
            .ok_or_else(|| anyhow::anyhow!("outgoing tangent missing"))?;
        let turn = angle_about(axis, incoming, outgoing);
        let setback = radius * (turn / 2.0).tan();
        let entry = sub(v, mul(incoming, setback));
        let exit = add(v, mul(outgoing, setback));
        let center = add(entry, mul(cross(axis, incoming), radius));
        let bend = TorusBend {
            center_mm: center,
            axis,
            major_radius_mm: radius,
            minor_radius_mm: 0.1,
            face_ids: vec![],
        };
        corners.push((entry, exit, bend, turn));
    }
    let mut segments = Vec::new();
    for i in 0..corners.len() {
        let (entry, exit, bend, turn) = &corners[(i + 1) % corners.len()];
        segments.push(
            line_evidence(corners[i].1, *entry).ok_or_else(|| anyhow::anyhow!("line missing"))?,
        );
        segments.push(
            bend_evidence(bend, *entry, *exit, *turn)
                .ok_or_else(|| anyhow::anyhow!("bend missing"))?,
        );
    }
    let authored = assemble(1, vec![], 0.1, corners[0].1, true, segments)
        .ok_or_else(|| anyhow::anyhow!("authored sweep missing"))?;
    let fragment = recover_round_sweep_fragment(&authored)?;
    let kernel = MonstertruckKernel::project_standard();
    let solid = kernel.evaluate(&fragment.model, fragment.root)?;
    let step = kernel.to_step(&solid)?;
    let source = ParsedStepSource::parse(step.as_bytes())?.exchange;
    let entities = &source.data[0].entities;
    let index = build_index(entities);
    let ids = crate::brep::manifold_solid_ids(entities);
    let [id] = ids.as_slice() else {
        anyhow::bail!("single solid expected")
    };
    let reduced = ReducedSolid::from_step(*id, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("source reduction unavailable"))?;
    let partition = Partition::from_solid(&reduced).ok_or_else(|| {
        anyhow::anyhow!(
            "source carriers not recognized: {:?}",
            reduced
                .faces()
                .map(|f| (&f.id, &f.surface))
                .collect::<Vec<_>>()
        )
    })?;
    Graph::from_partition(partition)
        .ok_or_else(|| anyhow::anyhow!("source section graph not certified"))?;
    let recovered = detect(*id, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("source graph not recognized"))?;
    assert!(recovered.closed);
    assert_eq!(recovered.segments.len(), 6);
    let rebuilt = recover_round_sweep_fragment(&recovered)?;
    let proof = SourceRecovery::new(*id, entities, &index, &kernel)
        .ok_or_else(|| anyhow::anyhow!("source certificate unavailable"))?;
    let candidate = kernel.evaluate(&rebuilt.model, rebuilt.root)?;
    let candidate_step = kernel.to_step(&candidate)?;
    let parsed = ParsedStepSource::parse(candidate_step.as_bytes())?.exchange;
    let candidate_entities = &parsed.data[0].entities;
    let candidate_ids = crate::brep::manifold_solid_ids(candidate_entities);
    let [candidate_id] = candidate_ids.as_slice() else {
        anyhow::bail!("candidate must have a single solid");
    };
    let comparison = crate::boundary_proof::compare_analytic_boundaries(
        *id,
        entities,
        *candidate_id,
        candidate_entities,
    );
    assert!(
        comparison.equivalent,
        "{comparison:#?}\nRecovered: {recovered:#?}"
    );
    assert!(proof.proves(&rebuilt.model, rebuilt.root, None)?);
    Ok(())
}
