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
fn interference_paths_are_invariant_under_permutation_and_reversal() {
    let p = |x, y| Point3::new(x, y, 0.0);
    let lines = vec![
        (p(0.0, 0.0), p(1.0, 0.0)),
        (p(1.0, 0.0), p(1.0, 1.0)),
        (p(1.0, 1.0), p(0.0, 1.0)),
        (p(0.0, 1.0), p(0.0, 0.0)),
        (p(2.0, 0.0), p(3.0, 0.0)),
        (p(3.0, 0.0), p(4.0, 0.0)),
    ];
    let baseline = construct_polylines(&lines);
    assert_eq!(baseline.len(), 2);
    assert_eq!(baseline.iter().map(|chain| chain.len()).sum::<usize>(), 8);
    for rotation in 0..lines.len() {
        let mut reordered = lines.clone();
        reordered.rotate_left(rotation);
        reordered.reverse();
        for segment in &mut reordered {
            std::mem::swap(&mut segment.0, &mut segment.1);
        }
        assert_eq!(construct_polylines(&reordered), baseline);
    }
}

#[test]
fn branch_junctions_are_not_arbitrarily_routed_into_one_chain() {
    let p = |x, y| Point3::new(x, y, 0.0);
    let a = p(-1.0, 0.0);
    let b = p(0.0, 0.0);
    let c = p(1.0, 0.0);
    let d = p(0.0, 1.0);
    let lines = vec![(b, c), (a, b), (b, d)];
    let baseline = construct_polylines(&lines);
    assert_eq!(baseline.len(), 3, "a 3-way junction must expose three arms");
    assert!(baseline.iter().all(|chain| chain.len() == 2));
    assert_eq!(
        baseline
            .iter()
            .filter(|chain| chain[0] == b || chain[1] == b)
            .count(),
        3
    );
    for permutation in [vec![(b, d), (a, b), (b, c)], vec![(c, b), (d, b), (b, a)]] {
        assert_eq!(construct_polylines(&permutation), baseline);
    }
}

#[test]
fn cycle_attached_to_branch_is_preserved_as_a_cycle_and_a_spur() {
    let p = |x, y| Point3::new(x, y, 0.0);
    let join = p(0.0, 0.0);
    let lines = vec![
        (join, p(1.0, 0.0)),
        (p(1.0, 0.0), p(1.0, 1.0)),
        (p(1.0, 1.0), p(0.0, 1.0)),
        (p(0.0, 1.0), join),
        (join, p(-1.0, 0.0)),
    ];
    let baseline = construct_polylines(&lines);
    assert_eq!(baseline.len(), 2);
    let cycle = baseline.iter().find(|chain| chain.len() == 5).unwrap();
    assert_eq!(cycle[0], join);
    assert_eq!(cycle[4], join);
    let spur = baseline.iter().find(|chain| chain.len() == 2).unwrap();
    assert!(spur.contains(&join));
    assert_eq!(
        construct_polylines(&lines.iter().rev().map(|&(a, b)| (b, a)).collect::<Vec<_>>()),
        baseline
    );
}

#[test]
fn quantized_node_representative_is_independent_of_input_order() {
    let p = |x| Point3::new(x, 0.0, 0.0);
    let a = p(0.0);
    let middle = p(1.0);
    let middle_near = p(1.0 + 0.25 * KERNEL_EPSILON);
    assert_eq!(
        PointIndex::from_point(middle),
        PointIndex::from_point(middle_near)
    );
    let c = p(2.0);
    let forward = construct_polylines(&[(a, middle), (middle_near, c)]);
    let backward = construct_polylines(&[(c, middle_near), (middle, a)]);
    assert_eq!(forward, backward);
    assert_eq!(forward.len(), 1);
    assert_eq!(forward[0].len(), 3);
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

#[test]
fn overlapping_mesh_strand_is_removed_before_ambiguous_endpoint_stitching() {
    let base = PolylineCurve(vec![
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    ]);
    let long = PolylineCurve(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]);
    let duplicated = PolylineCurve(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.25, 0.0, 0.0),
    ]);
    assert!(
        stitch_nearby_polylines(vec![base.clone(), long.clone(), duplicated.clone()], 1.0e-6)
            .is_none()
    );
    let clean = remove_redundant_interference_strands(vec![base, long, duplicated], 1.0e-6);
    assert_eq!(clean.len(), 2);
    let joined = stitch_nearby_polylines(clean, 1.0e-6)
        .expect("the genuine intersection path is no longer branched");
    assert_eq!(joined.len(), 1);
    assert_eq!(joined[0].len(), 3);
}

#[test]
fn distinct_meshed_branches_are_not_silently_deleted() {
    let base = PolylineCurve(vec![
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    ]);
    let branch0 = PolylineCurve(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]);
    let branch1 = PolylineCurve(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.25, 0.01, 0.0),
    ]);
    let clean = remove_redundant_interference_strands(vec![base, branch0, branch1], 1.0e-6);
    assert_eq!(clean.len(), 3);
    assert!(stitch_nearby_polylines(clean, 1.0e-6).is_none());
}
