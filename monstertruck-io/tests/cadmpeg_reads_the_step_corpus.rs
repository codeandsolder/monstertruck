//! The cadmpeg-backed STEP reader must read the same corpus ours does.
//!
//! `monstertruck-io` carries two STEP readers on purpose: our own on `step-p21`
//! (the default, and the measurement baseline the boolean kernel is
//! characterised against) and cadmpeg's behind the `cadmpeg` feature. Keeping
//! the second one compiled but unexercised would tell us nothing about whether
//! it could ever replace the first, which is the only reason it is here.
//!
//! So this row reads every in-repo STEP fixture through cadmpeg and asserts the
//! recovered *topology*, not merely that decoding returned `Ok`. A decoder that
//! parses a file and recovers no faces has not read it in any sense that matters
//! to a B-rep kernel -- and that was the real failure mode: on 2026-08-04
//! `occt-cylinder.step` decoded "successfully" while dropping its entire
//! `MANIFOLD_SOLID_BREP`, because one parameter-space `CIRCLE` failed
//! (cadmpeg/cadmpeg#79, fixed in cadmpeg/cadmpeg#83).
//!
//! The expected topology is read independently from the STEP entities in each
//! fixture. The cadmpeg decoder therefore never supplies its own oracle.

#![cfg(feature = "cadmpeg")]

use monstertruck_io::cadmpeg::step::decode_file_to_ir;
use std::path::{Path, PathBuf};

const FIXTURES: [&str; 10] = [
    "abc-0000.step",
    "abc-0006.step",
    "abc-0008.step",
    "abc-0035.step",
    "occt-assy.step",
    "occt-cone.step",
    "occt-cube.step",
    "occt-cylinder.step",
    "occt-sphere.step",
    "occt-torus.step",
];

fn occurrence_count(haystack: &[u8], needle: &[u8]) -> usize {
    haystack
        .windows(needle.len())
        .filter(|window| *window == needle)
        .count()
}

fn encoded_topology(path: &Path) -> (usize, usize, usize, usize) {
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    (
        occurrence_count(&bytes, b"MANIFOLD_SOLID_BREP("),
        occurrence_count(&bytes, b"ADVANCED_FACE("),
        occurrence_count(&bytes, b"EDGE_CURVE("),
        occurrence_count(&bytes, b"VERTEX_POINT("),
    )
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../resources/step")
        .join(name)
}

#[test]
fn cadmpeg_recovers_the_expected_topology_from_every_fixture() {
    let mut failures = Vec::new();

    for name in FIXTURES {
        let path = fixture(name);
        assert!(
            path.is_file(),
            "fixture {} is missing -- it moved, or a rename left this path behind",
            path.display(),
        );
        let expected = encoded_topology(&path);

        let ir = match decode_file_to_ir(&path) {
            Ok(ir) => ir,
            Err(error) => {
                failures.push(format!("{name}: decode failed: {error}"));
                continue;
            }
        };

        let got = (
            ir.model.bodies.len(),
            ir.model.faces.len(),
            ir.model.edges.len(),
            ir.model.vertices.len(),
        );
        if got != expected {
            failures.push(format!(
                "{name}: decoded bodies/faces/edges/vertices = {got:?}, but the STEP file encodes {expected:?}"
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "{} of {} fixtures did not match their encoded topology:\n  {}",
        failures.len(),
        FIXTURES.len(),
        failures.join("\n  "),
    );
}

/// The fixture that regressed, pinned on its own so a failure names the defect.
#[test]
fn the_cylinder_that_once_lost_its_whole_solid_now_keeps_it() {
    let ir = decode_file_to_ir(fixture("occt-cylinder.step")).expect("decode occt-cylinder");

    // Not "did it parse" -- did the MANIFOLD_SOLID_BREP survive. A cylinder is
    // three faces: the lateral surface and two caps.
    assert_eq!(
        ir.model.faces.len(),
        3,
        "the solid collapsed again: cadmpeg/cadmpeg#79 was one parameter-space \
         CIRCLE taking the entire MANIFOLD_SOLID_BREP with it"
    );
    assert_eq!(ir.model.bodies.len(), 1, "expected one body");
}
