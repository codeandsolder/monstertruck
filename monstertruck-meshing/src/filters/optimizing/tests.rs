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
fn out_of_range_bucket_coordinate_is_not_deduplicated() {
    let attrs = [Vector2::new(f64::MAX, 0.0), Vector2::new(f64::MAX, 0.0)];
    let blocks = create_blocks(&attrs, 1.0e-6);
    assert!(blocks.is_empty());
}
