//! Fast source-derived regression for exact SSI trim attachment. The public
//! difference and independent solid-angle occupancy oracle live in the
//! corresponding `boolean/ssi_radial_slot.rs` integration test.
use super::*;
use monstertruck_modeling::{Curve, Solid as ModelingSolid, Surface};

fn saved_operands() -> (ModelingSolid, ModelingSolid) {
    let host = serde_json::from_slice(include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/radial-slot-source/host.json"
    )))
    .expect("source-derived exact host B-rep");
    let cutter = serde_json::from_slice(include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/radial-slot-source/cutter.json"
    )))
    .expect("source-derived exact cutter B-rep");
    (host, cutter)
}

type RadialFaceSplits = (
    loops_store::InterfaceLoops<AltCurve<Curve, Surface>>,
    AltCurveShell<Curve, Surface>,
    AltCurveShell<Curve, Surface>,
);

fn split_original_faces(
    source: Shell<Point3, Curve, Surface>,
    cutter: Shell<Point3, Curve, Surface>,
    source_face: usize,
) -> RadialFaceSplits {
    let numeric =
        monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length;
    let mesh_seed = BooleanTolerance::new(numeric).mesh_seed;
    let polygon0 = source.triangulation(mesh_seed);
    let polygon1 = cutter.triangulation(mesh_seed);
    let geom0: AltCurveShell<Curve, Surface> = source.mapped(
        |point| *point,
        |curve| Alternative::FirstType(curve.clone()),
        Clone::clone,
    );
    let geom1: AltCurveShell<Curve, Surface> = cutter.mapped(
        |point| *point,
        |curve| Alternative::FirstType(curve.clone()),
        Clone::clone,
    );
    let loops =
        loops_store::create_loops_stores(&geom0, &polygon0, &geom1, &polygon1, numeric, true)
            .unwrap_or_else(|| panic!("source-derived radial slot face {source_face} SSI lost"));
    (loops, geom0, geom1)
}

#[test]
fn radial_slot_bore_ssi_endpoints_remain_on_original_exact_trims() {
    let (host, cutter) = saved_operands();
    // The first half of the narrow bore is exactly cut at its original axial
    // end trims. The box cutter extends 0.0001 mm farther in that direction.
    let bore = Shell::from(vec![host.boundaries()[0][30].clone()]);
    let original_bounds = bore[0]
        .vertex_iter()
        .map(|vertex| vertex.point().z)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), z| {
            (lo.min(z), hi.max(z))
        });
    let (loops, geom, _) = split_original_faces(bore, cutter.boundaries()[0].clone(), 30);
    let wires = &loops.geom[0][0];
    assert_eq!(
        wires.iter().map(|wire| wire.len()).collect::<Vec<_>>(),
        [4, 4]
    );
    for edge in wires.iter().flat_map(|wire| wire.iter()) {
        assert!(edge.is_geometric_consistent());
    }
    for vertex in wires.iter().flat_map(|wire| wire.vertex_iter()) {
        let z = vertex.point().z;
        assert!(z >= original_bounds.0 - 1.0e-6);
        assert!(z <= original_bounds.1 + 1.0e-6);
    }
    let numeric =
        monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length;
    let regions =
        divide_face::divide_faces(&geom, &loops.geom[0], loops.sampled[0].as_ref(), numeric)
            .expect("original bore remains divisible after exact SSI attachment");
    assert_eq!(regions.all_regions().len(), 2);
}

#[test]
fn radial_slot_curved_ssi_arcs_partition_the_flat_cutter_cap() {
    let (host, cutter) = saved_operands();
    // Both semicylinders cross the upper cutter cap, with intersection
    // leaders overrunning its side trims by sampled marching error.
    let source = Shell::from(vec![
        host.boundaries()[0][16].clone(),
        host.boundaries()[0][48].clone(),
    ]);
    let cap = Shell::from(vec![cutter.boundaries()[0][3].clone()]);
    let (loops, _, geom) = split_original_faces(source, cap, 16);
    assert_eq!(
        loops.geom[1][0]
            .iter()
            .map(|wire| wire.len())
            .collect::<Vec<_>>(),
        [5, 5]
    );
    let numeric =
        monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length;
    let regions =
        divide_face::divide_faces(&geom, &loops.geom[1], loops.sampled[1].as_ref(), numeric)
            .expect("curved intersection arcs must separate cap material regions");
    assert_eq!(regions.all_regions().len(), 2);
}
