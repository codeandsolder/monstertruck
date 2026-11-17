use super::common::solid_angle_contains;
use anyhow::Result;
use monstertruck_geometry::prelude::ParametricCurve;
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

    // Independently measured analytic STEP source, not a point inferred from
    // the Boolean result: the two slot-root planar/cylindrical SSI arcs have
    // radius 0.4826 mm about this exact source cylinder axis. Endpoints alone
    // are insufficient: the previous degree-1 spline export displaced the
    // interior by 0.0095 mm while retaining every topological endpoint.
    let source_axis_xy = (59.384_373_707_781_98_f64, -5.257_812_234_042_557_f64);
    let source_root_z = -2.590_787_765_957_25_f64;
    let source_radius = 0.4826_f64;
    let on_root_circle = |point: Point3| {
        (point.z - source_root_z).abs() <= 1.0e-6
            && ((point.x - source_axis_xy.0).hypot(point.y - source_axis_xy.1) - source_radius)
                .abs()
                <= 1.0e-6
    };
    let mut certified_root_edges = std::collections::HashSet::new();
    for edge in result.edge_iter() {
        let (front, back) = (edge.absolute_front().point(), edge.absolute_back().point());
        if !on_root_circle(front) || !on_root_circle(back) || front.distance(back) < 0.1 {
            continue;
        }
        if !certified_root_edges.insert(edge.id()) {
            continue;
        }
        let curve = edge.curve();
        let monstertruck_modeling::Curve::Intersection(intersection) = &curve else {
            anyhow::bail!("slot-root SSI arc lost its native two-carrier intersection identity");
        };
        let monstertruck_modeling::Curve::NurbsCurve(guide) = intersection.leader().as_ref() else {
            anyhow::bail!("slot-root SSI leader degraded to a piecewise-linear guide");
        };
        anyhow::ensure!(
            guide.degree() == 2,
            "slot-root SSI must retain its quadratic rational circle"
        );
        let (lo, hi) = curve.range_tuple();
        for sample in 0..=64 {
            let t = f64::mul_add(hi - lo, f64::from(sample) / 64.0, lo);
            let point = curve.evaluate(t);
            anyhow::ensure!(
                on_root_circle(point),
                "exact planar-cylinder SSI arc strayed from the source circle at sample {sample}:                  {point:?}, radius error={}",
                (point.x - source_axis_xy.0).hypot(point.y - source_axis_xy.1) - source_radius,
            );
        }
    }
    anyhow::ensure!(
        certified_root_edges.len() >= 2,
        "the Boolean result lost the two independent exact circular SSI root arcs"
    );
    Ok(())
}
