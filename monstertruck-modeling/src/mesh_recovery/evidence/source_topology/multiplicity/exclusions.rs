//! An embedded oriented shell cannot contain two coincident patches with
//! opposite winding. Exact simple planar polygon copies provide a witness
//! of this incompatibility; neither copy is removed or displaced.
use super::super::super::super::InputPolygon;
use super::equivalence::{polygon_key, simple_planar_polygon};
use crate::Point3;
use std::collections::BTreeMap;

pub(super) struct Exclusion {
    pub(super) classes: [usize; 2],
    pub(super) polygons: [usize; 2],
}

pub(super) fn opposed_faces(
    positions: &[Point3],
    polygons: &[InputPolygon],
    fragments: &[&[usize]],
    classes: &[usize],
) -> Vec<Exclusion> {
    // Keep a representative per geometric cycle and fragment class. A
    // class has the same oriented coverage for every instance, so an
    // overlap witness excludes every instance of the other class.
    let mut cycles = BTreeMap::<Vec<usize>, BTreeMap<usize, usize>>::new();
    for (&class, fragment) in classes.iter().zip(fragments) {
        for &id in *fragment {
            cycles
                .entry(polygon_key(&polygons[id]))
                .or_default()
                .entry(class)
                .or_insert(id);
        }
    }
    let mut pairs = BTreeMap::new();
    for (cycle, left) in &cycles {
        if cycle.len() < 3 {
            continue;
        }
        let mut reversed = cycle.clone();
        reversed[1..].reverse();
        if *cycle >= reversed {
            continue;
        }
        let Some(right) = cycles.get(&reversed) else {
            continue;
        };
        let Some(&representative) = left.values().next() else {
            continue;
        };
        if !simple_planar_polygon(positions, &polygons[representative]) {
            continue;
        }
        for (&a, &p) in left {
            for (&b, &q) in right {
                let mut classes = [a, b];
                classes.sort_unstable();
                let mut witness = [p, q];
                witness.sort_unstable();
                pairs.entry(classes).or_insert(witness);
            }
        }
    }
    pairs
        .into_iter()
        .map(|(classes, polygons)| Exclusion { classes, polygons })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::input_polygon_from_ids;

    #[test]
    fn exact_reversed_planar_cycles_exclude_classes_without_canceling()
    -> Result<(), Box<dyn std::error::Error>> {
        let positions = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 1e-10),
        ];
        let faces = [
            vec![0, 1, 2, 3],
            vec![2, 1, 0, 3],
            vec![1, 2, 3, 0],
            vec![0, 1, 2, 4],
            vec![0, 4, 2, 1],
        ];
        let polygons = faces
            .into_iter()
            .enumerate()
            .map(|(id, face)| input_polygon_from_ids(&positions, face, id, crate::KERNEL_EPSILON))
            .collect::<Result<Vec<_>, _>>()?;
        let result = opposed_faces(
            &positions,
            &polygons,
            &[&[0], &[1], &[2], &[3], &[4]],
            &[4, 5, 4, 6, 7],
        );
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].classes, [4, 5]);
        assert_eq!(result[0].polygons, [0, 1]);
        Ok(())
    }
}
