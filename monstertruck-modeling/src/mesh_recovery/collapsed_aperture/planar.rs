//! Constrained triangulation of planar domains with unchanged boundary IDs.
//! Reject intersecting, touching or incompletely recovered constraint loops.
use super::super::{edge_key, trace_oriented_boundary_cycles};
use robust::{Coord, orient2d};
use rustc_hash::FxHashMap;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

pub(super) fn turn(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let point = |[x, y]| Coord { x, y };
    orient2d(point(a), point(b), point(c))
}

pub(super) fn edges(ring: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    ring.iter()
        .copied()
        .zip(ring.iter().copied().cycle().skip(1))
}

pub(super) fn area_twice(points: &[[f64; 2]], ring: &[usize]) -> f64 {
    let p = points[ring[0]];
    edges(ring)
        .map(|(a, b)| turn(p, points[a], points[b]))
        .sum()
}

fn on_segment(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
    turn(a, b, p) == 0.0 && (0..2).all(|i| p[i] >= a[i].min(b[i]) && p[i] <= a[i].max(b[i]))
}

fn intersects(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let s = [turn(a, b, c), turn(a, b, d), turn(c, d, a), turn(c, d, b)];
    s.iter().any(|v| !v.is_finite())
        || ((s[0] > 0.0 && s[1] < 0.0 || s[0] < 0.0 && s[1] > 0.0)
            && (s[2] > 0.0 && s[3] < 0.0 || s[2] < 0.0 && s[3] > 0.0))
        || on_segment(a, b, c)
        || on_segment(a, b, d)
        || on_segment(c, d, a)
        || on_segment(c, d, b)
}

pub(super) fn simple(points: &[[f64; 2]], ring: &[usize]) -> bool {
    let sides = edges(ring).collect::<Vec<_>>();
    sides.len() >= 3
        && sides.iter().enumerate().all(|(i, &(a, b))| {
            let c = sides[(i + 1) % sides.len()].1;
            a != b
                && !on_segment(points[a], points[b], points[c])
                && !on_segment(points[b], points[c], points[a])
                && sides.iter().enumerate().skip(i + 2).all(|(j, &(c, d))| {
                    (i == 0 && j + 1 == sides.len())
                        || !intersects(points[a], points[b], points[c], points[d])
                })
        })
}

fn winding(points: &[[f64; 2]], ring: &[usize], p: [f64; 2]) -> Option<i32> {
    let mut count = 0_i32;
    for (a, b) in edges(ring) {
        let [a, b] = [points[a], points[b]];
        if on_segment(a, b, p) {
            return None;
        }
        let sign = turn(a, b, p);
        if a[1] <= p[1] && p[1] < b[1] && sign > 0.0 {
            count = count.checked_add(1)?;
        }
        if b[1] <= p[1] && p[1] < a[1] && sign < 0.0 {
            count = count.checked_sub(1)?;
        }
    }
    Some(count)
}

pub(super) fn inside_domain(points: &[[f64; 2]], loops: &[Vec<usize>], p: [f64; 2]) -> bool {
    loops
        .iter()
        .try_fold(0_i32, |n, ring| n.checked_add(winding(points, ring, p)?))
        == Some(1)
}

pub(super) fn disjoint(points: &[[f64; 2]], a: &[usize], b: &[usize]) -> bool {
    edges(a).all(|(i, j)| {
        edges(b).all(|(k, l)| !intersects(points[i], points[j], points[k], points[l]))
    })
}

pub(super) fn valid_domain(points: &[[f64; 2]], loops: &[Vec<usize>]) -> bool {
    let Some(outer) = loops.iter().find(|r| area_twice(points, r) > 0.0) else {
        return false;
    };
    loops.iter().all(|r| {
        simple(points, r)
            && (std::ptr::eq(r, outer)
                || (area_twice(points, r) < 0.0 && winding(points, outer, points[r[0]]) == Some(1)))
    }) && loops.iter().enumerate().all(|(i, r)| {
        loops.iter().skip(i + 1).all(|s| {
            disjoint(points, r, s)
                && (std::ptr::eq(r, outer)
                    || std::ptr::eq(s, outer)
                    || (winding(points, r, points[s[0]]) == Some(0)
                        && winding(points, s, points[r[0]]) == Some(0)))
        })
    })
}

// Oriented chain boundary with exactly two opposite uses for each
// interior edge. Vertices are not merged or simplified.
pub(super) fn boundary<'a>(polygons: impl Iterator<Item = &'a [usize]>) -> Option<Vec<Vec<usize>>> {
    let mut incidence = FxHashMap::<[usize; 2], (usize, i32)>::default();
    for polygon in polygons {
        for (a, b) in edges(polygon) {
            let (key, forward) = edge_key(a, b);
            let record = incidence.entry(key).or_default();
            record.0 += 1;
            record.1 += if forward { 1 } else { -1 };
        }
    }
    let mut next = FxHashMap::default();
    let mut incoming = FxHashMap::default();
    for ([a, b], (count, sign)) in incidence {
        let (a, b) = match (count, sign) {
            (2, 0) => continue,
            (1, 1) => (a, b),
            (1, -1) => (b, a),
            _ => return None,
        };
        if next.insert(a, b).is_some() || incoming.insert(b, a).is_some() {
            return None;
        }
    }
    if next.len() != incoming.len() || next.keys().any(|id| !incoming.contains_key(id)) {
        return None;
    }
    trace_oriented_boundary_cycles(next)
}

fn boundary_edges(loops: &[Vec<usize>]) -> Vec<(usize, usize)> {
    let mut result = loops
        .iter()
        .flat_map(|ring| edges(ring))
        .collect::<Vec<_>>();
    result.sort_unstable();
    result
}

pub(super) fn triangulate(points: &[[f64; 2]], loops: &[Vec<usize>]) -> Option<Vec<[usize; 3]>> {
    if !valid_domain(points, loops) {
        return None;
    }
    let mut ids = loops.iter().flatten().copied().collect::<Vec<_>>();
    ids.sort_unstable();
    ids.dedup();
    let local = ids
        .iter()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect::<FxHashMap<_, _>>();
    let vertices = ids
        .iter()
        .map(|&id| Point2::new(points[id][0], points[id][1]))
        .collect::<Vec<_>>();
    let constraints = loops
        .iter()
        .flat_map(|ring| edges(ring).map(|(a, b)| [local[&a], local[&b]]))
        .collect::<Vec<_>>();
    let count = constraints.len();
    let mut conflict = false;
    let cdt = ConstrainedDelaunayTriangulation::try_bulk_load_cdt(vertices, constraints, |_| {
        conflict = true
    })
    .ok()?;
    if conflict || cdt.num_constraints() != count || cdt.num_vertices() != ids.len() {
        return None;
    }
    let mut triangles = Vec::new();
    for face in cdt.inner_faces() {
        let mut tri = face.vertices().map(|v| ids[v.index()]);
        let [a, b, c] = tri.map(|id| points[id]);
        let centroid = [
            a[0] + ((b[0] - a[0]) + (c[0] - a[0])) / 3.0,
            a[1] + ((b[1] - a[1]) + (c[1] - a[1])) / 3.0,
        ];
        if inside_domain(points, loops, centroid) {
            if turn(a, b, c) < 0.0 {
                tri.swap(1, 2);
            }
            triangles.push(tri);
        }
    }
    let recovered = boundary(triangles.iter().map(<[usize; 3]>::as_slice))?;
    (boundary_edges(&recovered) == boundary_edges(loops)).then_some(triangles)
}
