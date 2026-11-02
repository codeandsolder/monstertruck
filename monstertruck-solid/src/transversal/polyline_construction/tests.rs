use super::*;

#[test]
fn construct_polylines_positive0() {
    let lines = vec![
        (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 0.0, 1.0)),
        (Point3::new(0.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)),
        (Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)),
        (Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(1.0, 1.0, 1.0), Point3::new(0.0, 0.0, 0.0)),
        (Point3::new(1.0, 0.0, 1.0), Point3::new(0.0, 1.0, 1.0)),
    ];
    let polyline = construct_polylines(&lines);
    assert_eq!(polyline.len(), 1);
    assert_eq!(polyline[0].len(), 9);

    let mut sign = None;
    for line in polyline[0].windows(2) {
        let a = line[0][0] + line[0][1] * 2.0 + line[0][2] * 4.0;
        let b = line[1][0] + line[1][1] * 2.0 + line[1][2] * 4.0;
        let x = b - a;
        assert!(f64::abs(x) == 1.0 || f64::abs(x) == 7.0);
        let s = f64::signum(x * (x - 2.0) * (x + 2.0));
        if let Some(sign) = sign {
            assert!(s == sign);
        } else {
            sign = Some(s);
        }
    }
}

#[test]
fn construct_polylines_positive1() {
    let lines = vec![
        (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(1.0, 0.0, 1.0), Point3::new(1.0, 1.0, 1.0)),
        (Point3::new(1.0, 1.0, 1.0), Point3::new(0.0, 1.0, 1.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
        (Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)),
        (Point3::new(0.0, 1.0, 1.0), Point3::new(0.0, 0.0, 1.0)),
    ];
    let polyline = construct_polylines(&lines);
    assert_eq!(polyline.len(), 2);
    assert_eq!(polyline[0].len(), 5);
    assert_eq!(polyline[1].len(), 5);
}

#[test]
fn construct_polylines_positive2() {
    let lines = vec![
        (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 0.0, 1.0)),
        (Point3::new(0.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)),
        (Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(1.0, 1.0, 1.0), Point3::new(0.0, 0.0, 0.0)),
        (Point3::new(1.0, 0.0, 1.0), Point3::new(0.0, 1.0, 1.0)),
    ];
    let polyline = construct_polylines(&lines);
    assert_eq!(polyline.len(), 1);
    assert_eq!(polyline[0].len(), 9);

    let mut sign = None;
    for line in polyline[0].windows(2) {
        let a = line[0][0] + line[0][1] * 2.0 + line[0][2] * 4.0;
        let b = line[1][0] + line[1][1] * 2.0 + line[1][2] * 4.0;
        let x = b - a;
        assert!(f64::abs(x) == 1.0 || f64::abs(x) == 7.0);
        let s = f64::signum(x * (x - 2.0) * (x + 2.0));
        if let Some(sign) = sign {
            assert!(s == sign);
        } else {
            sign = Some(s);
        }
    }
}

#[test]
fn construct_polylines_positive3() {
    let lines = vec![
        (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 0.0, 1.0)),
        (Point3::new(0.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)),
        (Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(1.0, 0.0, 1.0), Point3::new(0.0, 1.0, 1.0)),
    ];
    let polyline = construct_polylines(&lines);
    assert_eq!(polyline.len(), 1);
    assert_eq!(polyline[0].len(), 8);

    let mut sign = None;
    for line in polyline[0].windows(2) {
        let a = line[0][0] + line[0][1] * 2.0 + line[0][2] * 4.0;
        let b = line[1][0] + line[1][1] * 2.0 + line[1][2] * 4.0;
        let x = b - a;
        assert!(f64::abs(x) == 1.0);
        let s = f64::signum(x * (x - 2.0) * (x + 2.0));
        if let Some(sign) = sign {
            assert!(s == sign);
        } else {
            sign = Some(s);
        }
    }
}

#[test]
fn duplicate_interference_segments_are_set_like() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);

    for lines in [vec![(a, b), (a, b)], vec![(a, b), (b, a)]] {
        let polylines = construct_polylines(&lines);
        assert_eq!(polylines.len(), 1);
        assert_eq!(polylines[0].len(), 2);
        let front = polylines[0][0];
        let back = polylines[0][polylines[0].len() - 1];
        assert!(
            (front.near(&a) && back.near(&b)) || (front.near(&b) && back.near(&a)),
            "deduplicated segment must preserve the two geometric endpoints"
        );
    }
}

#[test]
fn duplicate_segment_does_not_turn_open_chain_into_closed_walk() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let c = Point3::new(2.0, 0.0, 0.0);
    let lines = vec![(a, b), (b, c), (c, b)];

    let polylines = construct_polylines(&lines);
    assert_eq!(polylines.len(), 1);
    assert_eq!(polylines[0].len(), 3);
    assert!(!polylines[0][0].near(&polylines[0][polylines[0].len() - 1]));
}

#[test]
fn nearby_endpoint_fragments_stitch_with_caller_tolerance() {
    let gap = 2.25e-6;
    let chains = vec![
        PolylineCurve(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]),
        PolylineCurve(vec![
            Point3::new(1.0 + gap, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ]),
        PolylineCurve(vec![
            Point3::new(2.0 + gap, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ]),
    ];
    let stitched = stitch_nearby_polylines(chains, 4.0e-6).expect("unambiguous endpoint stitching");
    assert_eq!(stitched.len(), 1);
    assert_eq!(stitched[0].len(), 4);
}

#[test]
fn nearby_endpoint_stitching_rejects_ambiguous_branches() {
    let chains = vec![
        PolylineCurve(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]),
        PolylineCurve(vec![
            Point3::new(1.0 + 1.0e-6, 1.0e-6, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ]),
        PolylineCurve(vec![
            Point3::new(1.0 + 1.0e-6, -1.0e-6, 0.0),
            Point3::new(2.0, -1.0, 0.0),
        ]),
    ];
    assert!(stitch_nearby_polylines(chains, 4.0e-6).is_none());
}

#[test]
fn non_finite_interference_segments_are_ignored() {
    let finite = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    for bad in [
        Point3::new(f64::NAN, 0.0, 0.0),
        Point3::new(f64::INFINITY, 0.0, 0.0),
        Point3::new(0.0, f64::NEG_INFINITY, 0.0),
    ] {
        let polylines = construct_polylines(&[finite, (finite.1, bad)]);
        assert_eq!(polylines.len(), 1);
        assert_eq!(polylines[0].len(), 2);
    }
}

#[test]
fn ordinary_point_index_hash_matches_legacy_bucket_hash() {
    use rustc_hash::FxHasher;
    use std::hash::{Hash, Hasher};

    let index = PointIndex::from_point(Point3::new(1.25, -2.5, 0.125))
        .expect("finite ordinary point must be indexable");
    let buckets = match index.0 {
        [
            CoordinateIndex::Bucket(x),
            CoordinateIndex::Bucket(y),
            CoordinateIndex::Bucket(z),
        ] => [x, y, z],
        _ => panic!("ordinary point unexpectedly used exact-coordinate fallback"),
    };

    let mut current = FxHasher::default();
    index.hash(&mut current);
    let mut legacy = FxHasher::default();
    buckets.hash(&mut legacy);
    assert_eq!(current.finish(), legacy.finish());
}

#[test]
fn large_finite_interference_coordinates_are_preserved_and_stitched() {
    let x = f64::MAX;
    let a = Point3::new(x, 0.0, 0.0);
    let b = Point3::new(x, 1.0, 0.0);
    let b_near = Point3::new(x, 1.0 + 0.5 * KERNEL_EPSILON, 0.0);
    let c = Point3::new(x, 2.0, 0.0);

    let polylines = construct_polylines(&[(a, b), (b_near, c)]);
    assert_eq!(polylines.len(), 1);
    assert_eq!(polylines[0].len(), 3);

    let front = polylines[0][0];
    let back = polylines[0][2];
    assert!(
        (front.near(&a) && back.near(&c)) || (front.near(&c) && back.near(&a)),
        "large finite coordinates must not make interference segments disappear"
    );
}
