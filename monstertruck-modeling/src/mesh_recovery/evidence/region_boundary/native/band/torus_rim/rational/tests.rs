use super::*;
use crate::{BsplineCurve, KnotVector, Vector4};

fn frame() -> Frame {
    Frame {
        center: Point3::new(0.0, 0.0, 0.0),
        basis: [Vector3::unit_x(), Vector3::unit_y(), Vector3::unit_z()],
        major: 3.0,
        varying: 0,
        direction: 1.0,
        side: 1.0,
        resolution: 1e-8,
    }
}

fn torus([u, v]: [f64; 2]) -> [f64; 3] {
    let (a, b) = u.sin_cos();
    let (c, d) = v.sin_cos();
    [(3.0 + d) * b, (3.0 + d) * a, c]
}
fn control(p: [f64; 3]) -> Control { [p[0], p[1], p[2], 1.0].map(Interval::point) }

fn polygonal_rim(points: &[[f64; 2]]) -> Option<Vec<Piece>> {
    let f = frame();
    let mut result = Vec::new();
    let mut reference = points.first()?.map(Interval::point);
    for (i, pair) in points.windows(2).enumerate() {
        let values = pair.iter().map(|&p| control(torus(p))).collect();
        reference = append(
            f,
            &Controls {
                values,
                parameter: [0.0, 1.0],
            },
            reference,
            i,
            MAX_DEPTH,
            &mut result,
        )?;
    }
    Some(result)
}

#[test]
fn actual_curve_cycle_proves_backtracking_and_rejects_crossing() {
    let simple = [
        [0.0, 1.0],
        [0.8, 1.0],
        [0.65, 1.35],
        [1.6, 1.35],
        [1.45, 1.0],
        [2.5, 1.0],
        [3.5, 1.0],
        [4.5, 1.0],
        [5.5, 1.0],
        [TAU, 1.0],
    ];
    let pieces = polygonal_rim(&simple).expect("regular chart pieces");
    assert!(simplicity::cycle(frame(), &pieces, true).is_ok());
    // Proper crossing in the lifted chart, with the same total winding.
    let crossed = [
        [0.0, 1.0],
        [0.4, 1.4],
        [0.05, 1.4],
        [0.45, 0.9],
        [2.5, 1.0],
        [3.5, 1.0],
        [4.5, 1.0],
        [5.5, 1.0],
        [TAU, 1.0],
    ];
    let pieces = polygonal_rim(&crossed).expect("regular chart pieces");
    assert!(simplicity::cycle(frame(), &pieces, true).is_err());
}

#[test]
fn contractible_crossing_reports_original_edge_pair() {
    let crossed = [[0.2, 0.2], [0.8, 0.8], [0.2, 0.8], [0.8, 0.2], [0.2, 0.2]];
    let pieces = polygonal_rim(&crossed).expect("regular chart pieces");
    let pairs =
        simplicity::cycle(frame(), &pieces, false).expect_err("non-simple contractible loop");
    assert!(pairs.contains(&[0, 2]));
}

#[test]
fn interval_extraction_covers_nonclamped_knots_and_tiny_active_spans()
-> Result<(), Box<dyn std::error::Error>> {
    let curve = BsplineCurve::new(
        KnotVector::from(vec![
            0.0,
            0.0,
            0.0,
            0.25,
            0.250_000_000_01,
            0.75,
            1.0,
            1.0,
            1.0,
        ]),
        vec![
            Vector4::new(3.0, 0.0, 1.0, 1.0),
            Vector4::new(2.0, 1.0, 0.0, 0.5),
            Vector4::new(1.0, 2.0, 1.0, 1.0),
            Vector4::new(0.0, 3.0, 0.0, 0.75),
            Vector4::new(-1.0, 4.0, 1.0, 1.0),
            Vector4::new(-2.0, 4.0, 0.0, 1.0),
        ],
    )?;
    let certified = controls::spans(&curve).ok_or("interval spans")?;
    let ordinary = monstertruck_geometry::proof_bezier_spans(&curve).ok_or("ordinary spans")?;
    assert_eq!(certified.len(), 4);
    for (bounds, approximate) in certified.iter().zip(ordinary) {
        for (enclosure, coefficient) in bounds.values.iter().zip(approximate.controls) {
            for (bound, value) in enclosure.iter().zip(coefficient) {
                assert!(
                    bound.contains(value) || bound.sub(Interval::point(value)).magnitude() < 1e-12
                );
            }
        }
    }
    Ok(())
}
