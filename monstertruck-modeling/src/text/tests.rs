//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;
use crate::Vector3;

const fn placement(scale: f64, y_flip: bool, z: f64) -> Placement {
    Placement {
        scale,
        y_flip,
        z,
        closure_tolerance: 1e-7,
        offset: [0.0; 2],
    }
}

/// Minimal synthetic font test: builds contours manually via the collector.
#[test]
fn contour_collector_line_segments() {
    let mut c = ContourCollector::new();
    // Square contour.
    c.move_to(0.0, 0.0);
    c.line_to(100.0, 0.0);
    c.line_to(100.0, 100.0);
    c.line_to(0.0, 100.0);
    c.close();
    assert_eq!(c.contours.len(), 1);
    // 3 explicit line_to segments; the close back to start is handled
    // implicitly by `contour_to_wire` when assembling the wire.
    assert_eq!(c.contours[0].2.len(), 3);
}

#[test]
fn contour_to_wire_square() {
    let mut c = ContourCollector::new();
    c.move_to(0.0, 0.0);
    c.line_to(1.0, 0.0);
    c.line_to(1.0, 1.0);
    c.line_to(0.0, 1.0);
    c.close();

    let (sx, sy, segs) = &c.contours[0];
    let wire = contour_to_wire(*sx, *sy, segs, placement(1.0, false, 0.0)).unwrap();
    assert!(wire.is_closed());
    assert_eq!(wire.len(), 4);
}

#[test]
fn contour_to_wire_with_bezier() {
    let mut c = ContourCollector::new();
    c.move_to(0.0, 0.0);
    c.quad_to(0.5, 1.0, 1.0, 0.0);
    c.line_to(0.0, 0.0);
    c.close();

    let (sx, sy, segs) = &c.contours[0];
    let wire = contour_to_wire(*sx, *sy, segs, placement(1.0, false, 0.0)).unwrap();
    assert!(wire.is_closed());
    assert_eq!(wire.len(), 2);
}

#[test]
fn contour_to_wire_cubic() {
    let mut c = ContourCollector::new();
    c.move_to(0.0, 0.0);
    c.curve_to(0.3, 1.0, 0.7, 1.0, 1.0, 0.0);
    c.line_to(0.0, 0.0);
    c.close();

    let (sx, sy, segs) = &c.contours[0];
    let wire = contour_to_wire(*sx, *sy, segs, placement(1.0, false, 0.0)).unwrap();
    assert!(wire.is_closed());
    assert_eq!(wire.len(), 2);
}

#[test]
fn contour_y_flip() {
    let pt = placement(0.01, true, 5.0).point(10.0, 20.0);
    assert!((pt.x - 0.1).abs() < 1e-10);
    assert!((pt.y - (-0.2)).abs() < 1e-10);
    assert!((pt.z - 5.0).abs() < 1e-10);
}

#[test]
fn contour_no_y_flip() {
    let pt = placement(0.01, false, 0.0).point(10.0, 20.0);
    assert!((pt.x - 0.1).abs() < 1e-10);
    assert!((pt.y - 0.2).abs() < 1e-10);
}

#[test]
fn multiple_contours() {
    let mut c = ContourCollector::new();
    // Outer square.
    c.move_to(0.0, 0.0);
    c.line_to(10.0, 0.0);
    c.line_to(10.0, 10.0);
    c.line_to(0.0, 10.0);
    c.close();
    // Inner square (hole).
    c.move_to(2.0, 2.0);
    c.line_to(8.0, 2.0);
    c.line_to(8.0, 8.0);
    c.line_to(2.0, 8.0);
    c.close();

    assert_eq!(c.contours.len(), 2);

    let wires: Vec<Wire> = c
        .contours
        .iter()
        .map(|(sx, sy, segs)| contour_to_wire(*sx, *sy, segs, placement(0.1, false, 0.0)).unwrap())
        .collect();

    assert_eq!(wires.len(), 2);
    assert!(wires.iter().all(|w| w.is_closed()));
}

#[test]
fn degenerate_contour_filtered() {
    assert!(is_degenerate_contour(
        5.0,
        5.0,
        &[Segment::Line(5.0, 5.0), Segment::Line(5.0, 5.0)]
    ));
    assert!(!is_degenerate_contour(0.0, 0.0, &[Segment::Line(1.0, 0.0)]));
}

fn font(data: &'static [u8]) -> Font { Font::new(data, 0).expect("valid synthetic font") }

const QUADRATIC: &[u8] = include_bytes!("fixtures/quadratic.ttf");
const CUBIC: &[u8] = include_bytes!("fixtures/cubic.otf");
const COLLECTION: &[u8] = include_bytes!("fixtures/collection.ttc");
const VARIABLE: &[u8] = include_bytes!("fixtures/variable.ttf");

fn glyph(font: &Font, ch: char) -> GlyphId {
    font.charmap().map_unicode(ch).expect("fixture character")
}

fn check_translated(original: &Wire, translated: &Wire, offset: Vector3) {
    assert_eq!(original.len(), translated.len());
    for (a, b) in original.edge_iter().zip(translated.edge_iter()) {
        // Check the interior of every curve, not merely its vertices:
        // control points must undergo the same affine transformation.
        for t in [0.0, 0.125, 0.5, 0.875, 1.0] {
            let expected = a.curve().evaluate(t) + offset;
            assert!(b.curve().evaluate(t).distance(expected) < 1e-10);
        }
    }
}

#[test]
fn opentype_quadratic_and_cubic_outlines_preserve_curve_degree() {
    for (data, degree) in [(QUADRATIC, 2), (CUBIC, 3)] {
        let font = font(data);
        let wires = glyph_profile(&font, glyph(&font, 'A'), &TextOptions::default()).unwrap();
        assert_eq!(wires.len(), 1);
        assert!(wires[0].is_closed());
        assert!(wires[0].edge_iter().any(|edge| {
            matches!(edge.curve(), Curve::BsplineCurve(curve) if curve.degree() == degree)
        }));
        let hole = glyph_profile(&font, glyph(&font, 'O'), &TextOptions::default()).unwrap();
        assert_eq!(hole.len(), 2);
        assert!(hole.iter().all(Wire::is_closed));
    }
}

#[test]
fn text_advances_translate_all_curve_geometry() {
    for data in [QUADRATIC, CUBIC] {
        let font = font(data);
        let opts = TextOptions {
            scale: Some(0.002),
            y_flip: false,
            z: 7.0,
            ..Default::default()
        };
        let single = glyph_profile(&font, glyph(&font, 'A'), &opts).unwrap();
        let repeated = text_profile(&font, "A A", &opts).unwrap();
        assert_eq!(repeated.len(), 2);
        check_translated(&single[0], &repeated[0], Vector3::new(0.0, 0.0, 0.0));
        check_translated(&single[0], &repeated[1], Vector3::new(1.7, 0.0, 0.0));
    }
}

#[test]
fn positioned_glyphs_apply_both_offsets_before_y_flip() {
    let font = font(QUADRATIC);
    let id = glyph(&font, 'A');
    let opts = TextOptions::default();
    let single = glyph_profile(&font, id, &opts).unwrap();
    let positioned = positioned_glyph_profiles(
        &font,
        [PositionedGlyph {
            glyph_id: id,
            position: [-130.0, 275.0],
        }],
        &opts,
    )
    .unwrap();
    assert_eq!(positioned.len(), 1);
    check_translated(&single[0], &positioned[0], Vector3::new(-0.13, -0.275, 0.0));
}

#[test]
fn collection_and_variable_instances_preserve_advance_and_geometry() {
    let regular = Font::new(COLLECTION, 0).unwrap();
    let bold = Font::new(COLLECTION, 1).unwrap();
    assert!(Font::new(COLLECTION, 2).is_none());
    let variable = font(VARIABLE);
    let heavy = variable
        .instance_builder()
        .variations([("wght", 900.0)])
        .build();
    for (font, expected_width, expected_advance) in [
        (regular, 0.5, 0.6),
        (bold, 0.75, 0.9),
        (variable, 0.5, 0.6),
        (heavy, 0.75, 0.9),
    ] {
        let opts = TextOptions {
            y_flip: false,
            ..Default::default()
        };
        let wires = text_profile(&font, "AA", &opts).unwrap();
        assert_eq!(wires.len(), 2);
        let width = wires[0]
            .edge_iter()
            .map(|edge| edge.front().point().x)
            .fold(0.0_f64, f64::max);
        assert!((width - expected_width).abs() < 1e-10);
        check_translated(
            &wires[0],
            &wires[1],
            Vector3::new(expected_advance, 0.0, 0.0),
        );
    }
}

#[test]
fn empty_and_invalid_placement_are_distinguished() {
    let font = font(QUADRATIC);
    assert!(glyph_profile(&font, glyph(&font, ' '), &TextOptions::default()).is_err());
    assert!(
        text_profile(&font, "   ", &TextOptions::default())
            .unwrap()
            .is_empty()
    );
    for scale in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let opts = TextOptions {
            scale: Some(scale),
            ..Default::default()
        };
        assert_eq!(
            text_profile(&font, "A", &opts).unwrap_err(),
            Error::InvalidTextOptions
        );
    }
    assert!(!is_degenerate_contour(
        0.0,
        0.0,
        &[Segment::Quad(2.0, 3.0, 0.0, 0.0)]
    ));
}

#[test]
fn closed_bezier_and_repeated_endpoint_keep_valid_topology() {
    // The cubic is nonconstant even though both endpoints and its
    // midpoint coincide; exact subdivision must find another interior point.
    let segments = [
        Segment::Cubic(1.0, 2.0, -1.0, -2.0, 0.0, 0.0),
        Segment::Line(0.0, 0.0),
    ];
    let wire = contour_to_wire(0.0, 0.0, &segments, placement(1.0, false, 0.0)).unwrap();
    assert!(wire.is_closed());
    assert_eq!(wire.len(), 2);
    assert!(wire.edge_iter().all(|edge| edge.front() != edge.back()));
    let segments = [
        Segment::Line(1.0, 0.0),
        Segment::Line(1.0, 1.0),
        Segment::Line(0.0, 0.0),
        Segment::Line(0.0, 0.0),
    ];
    let wire = contour_to_wire(0.0, 0.0, &segments, placement(1.0, false, 0.0)).unwrap();
    assert_eq!(wire.len(), 3);
    assert!(wire.is_closed());
}
