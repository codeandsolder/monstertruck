use super::common::solid_angle_contains;
use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::Solid;
use monstertruck_solid::{GeometryPolicy, difference};

#[test]
fn source_4a6356_upper_round_exact_interface_closes() -> Result<()> {
    let source: Solid =
        serde_json::from_slice(include_bytes!("../fixtures/ssi-4a6356/stage4-solid.json"))?;
    let cutter: Solid = serde_json::from_slice(include_bytes!(
        "../fixtures/ssi-4a6356/upper-cutter-solid.json"
    ))?;
    let result = difference(&source, &cutter, GeometryPolicy::standard(1.0e-5))?;
    assert_eq!(result.boundaries().len(), 1);
    assert_eq!(result.boundaries()[0].len(), 47);
    assert!(result.boundaries()[0].check_solid_boundary().is_ok());
    assert!(result.is_geometric_consistent());
    // An independent solid-angle oracle (rather than the production ray
    // classifier) must reproduce A \ B in source, cutter and union bounds.
    // Tessellate each input before comparing to prevent shared mutable state.
    let source_mesh = source.triangulation(0.005).to_polygon();
    let cutter_mesh = cutter.triangulation(0.005).to_polygon();
    let result_mesh = result.triangulation(0.005).to_polygon();
    let volumes = [
        source_mesh.volume().abs(),
        cutter_mesh.volume().abs(),
        result_mesh.volume().abs(),
    ];
    let slack = (volumes[0] + volumes[1] + 1.0) * 0.001;
    anyhow::ensure!(
        volumes[2] <= volumes[0] + slack,
        "difference exceeds source volume: {volumes:?}"
    );
    anyhow::ensure!(
        volumes[2] + slack >= (volumes[0] - volumes[1]).max(0.0),
        "difference loses more than cutter volume: {volumes:?}"
    );
    let source_box = source_mesh.bounding_box();
    let cutter_box = cutter_mesh.bounding_box();
    let mut both = source_box;
    both += cutter_box;
    let boxes = [source_box, cutter_box, both];
    let fractions = [0.137, 0.311, 0.509, 0.707, 0.887];
    let mut checked = 0;
    let mut inside_source = 0;
    let mut inside_cutter = 0;
    let mut overlap = 0;
    for bounds in boxes {
        let min = bounds.min();
        let span = bounds.diagonal();
        for fx in fractions {
            for fy in fractions {
                for fz in fractions {
                    let point = Point3::new(
                        span.x.mul_add(fx, min.x),
                        span.y.mul_add(fy, min.y),
                        span.z.mul_add(fz, min.z),
                    );
                    let in_source = solid_angle_contains(&source_mesh, point);
                    let in_cutter = solid_angle_contains(&cutter_mesh, point);
                    let actual = solid_angle_contains(&result_mesh, point);
                    anyhow::ensure!(
                        actual == (in_source && !in_cutter),
                        "difference occupancy mismatch at {point:?}: source={in_source}, cutter={in_cutter}, result={actual}"
                    );
                    checked += 1;
                    inside_source += usize::from(in_source && !in_cutter);
                    inside_cutter += usize::from(!in_source && in_cutter);
                    overlap += usize::from(in_source && in_cutter);
                }
            }
        }
    }
    anyhow::ensure!(
        checked == 375 && inside_source > 0 && overlap > 0,
        "insufficient independent occupancy coverage: checked={checked} source_only={inside_source} cutter_only={inside_cutter} overlap={overlap}"
    );
    Ok(())
}
