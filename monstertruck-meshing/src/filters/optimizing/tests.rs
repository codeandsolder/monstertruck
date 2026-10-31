//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;

fn into_vertices(iter: &[usize]) -> Vec<Vertex> { iter.iter().map(|i| i.into()).collect() }

#[test]
fn degenerate_polygon_test() {
    let poly = into_vertices(&[0, 1, 2, 0, 3, 4, 5, 6, 3, 7, 8, 9]);
    let polys = split_into_nondegenerate(poly);
    assert_eq!(polys[0], into_vertices(&[0, 1, 2]));
    assert_eq!(polys[1], into_vertices(&[3, 4, 5, 6]));
    assert_eq!(polys[2], into_vertices(&[3, 7, 8, 9, 0]));
}

#[test]
fn unbucketable_coordinates_keep_deduplication_semantics() {
    let min_bucket = -(2.0_f64).powi(63);
    let attrs = [
        Vector2::new(0.0, 0.0),
        Vector2::new(10.0, 0.0),
        Vector2::new(f64::MAX, 1.0),
        Vector2::new(f64::MAX, 1.0),
        Vector2::new(min_bucket, 2.0),
        Vector2::new(min_bucket, 2.0),
    ];

    let representatives = sub_put_together_same_attrs(&attrs, 1.0);

    assert_eq!(representatives[0], 0);
    assert_eq!(representatives[1], 1);
    assert_eq!(representatives[2], representatives[3]);
    assert!([2, 3].contains(&representatives[2]));
    assert_eq!(representatives[4], representatives[5]);
    assert!([4, 5].contains(&representatives[4]));
}
