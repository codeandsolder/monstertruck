use super::common::solid_angle_contains;
use anyhow::Result;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::Solid;
use monstertruck_solid::{GeometryPolicy, difference};

/// Source-derived slotted revolution: the cutter extends past one exact trim,
/// while two curved faces must split the opposite cutter cap. Watertightness
/// alone is not enough: check material occupancy with an independent oracle.
#[test]
fn source_radial_slot_exact_difference_preserves_material_and_manifold() -> Result<()> {
    let source: Solid =
        serde_json::from_slice(include_bytes!("../fixtures/radial-slot-source/host.json"))?;
    let cutter: Solid =
        serde_json::from_slice(include_bytes!("../fixtures/radial-slot-source/cutter.json"))?;
    let result = difference(&source, &cutter, GeometryPolicy::standard(1.0e-5))?;
    assert_eq!(result.boundaries().len(), 1);
    assert!(result.boundaries()[0].check_solid_boundary().is_ok());
    assert!(result.is_geometric_consistent());

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
        volumes[2] <= volumes[0] + slack
            && volumes[2] + slack >= (volumes[0] - volumes[1]).max(0.0),
        "difference violated universal volume bounds: {volumes:?}"
    );

    let source_box = source_mesh.bounding_box();
    let cutter_box = cutter_mesh.bounding_box();
    let mut union_box = source_box;
    union_box += cutter_box;
    let fractions = [0.137, 0.311, 0.509, 0.707, 0.887];
    let mut checked = 0;
    let mut surviving = 0;
    let mut removed = 0;
    for bounds in [source_box, cutter_box, union_box] {
        let min = bounds.min();
        let span = bounds.diagonal();
        for x in fractions {
            for y in fractions {
                for z in fractions {
                    let point = Point3::new(
                        span.x.mul_add(x, min.x),
                        span.y.mul_add(y, min.y),
                        span.z.mul_add(z, min.z),
                    );
                    let inside_source = solid_angle_contains(&source_mesh, point);
                    let inside_cutter = solid_angle_contains(&cutter_mesh, point);
                    let actual = solid_angle_contains(&result_mesh, point);
                    anyhow::ensure!(
                        actual == (inside_source && !inside_cutter),
                        "incorrect radial-slot material at {point:?}: source={inside_source}, cutter={inside_cutter}, result={actual}"
                    );
                    checked += 1;
                    surviving += usize::from(inside_source && !inside_cutter);
                    removed += usize::from(inside_source && inside_cutter);
                }
            }
        }
    }
    anyhow::ensure!(
        checked == 375 && surviving > 0 && removed > 0,
        "radial-slot occupancy sample missed material or removed volume: surviving={surviving}, removed={removed}"
    );
    Ok(())
}
