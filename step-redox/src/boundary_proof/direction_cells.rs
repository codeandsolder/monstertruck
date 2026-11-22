//! A quantization grid is a representation, not a geometric discontinuity.
//! Two plane normals on opposite sides of a cell boundary can span less
//! than one direction quantum, precisely the width of an ordinary key cell.
//! Match such carriers only with unique correspondence, unchanged offset
//! and sense keys, and an identical complete remapped boundary signature.
use super::{
    AnalyticBoundarySignature, BoundaryChainKey, DIRECTION_Q, OrientedSupportKey,
    boundary_surface_support, canonical_boundary_direction, oriented_support_key,
    strict_closed_shell_faces,
};
use crate::brep::{self, SurfaceSupport};
use crate::math3::normalize;
use crate::step_graph::build_index;
use std::collections::{BTreeMap, BTreeSet};
use step_p21::ast::EntityInstance;

type Normals = BTreeMap<OrientedSupportKey, Vec<[f64; 3]>>;

pub(super) fn plane_normals(solid: u64, entities: &[EntityInstance]) -> Option<Normals> {
    let index = build_index(entities);
    let mut normals = Normals::new();
    for face_id in strict_closed_shell_faces(solid, entities, &index).ok()? {
        let key = oriented_support_key(face_id, entities, &index)?;
        if matches!(key, OrientedSupportKey::Plane { .. }) {
            let surface = brep::face_surface(face_id, entities, &index)?;
            let SurfaceSupport::Plane(plane) = boundary_surface_support(surface, entities, &index)
            else {
                return None;
            };
            let (normal, _) = canonical_boundary_direction(normalize(plane.normal, 0.0)?)?;
            normals.entry(key).or_default().push(normal);
        }
    }
    Some(normals)
}

fn adjacent_plane_keys(a: &OrientedSupportKey, b: &OrientedSupportKey) -> bool {
    match (a, b) {
        (
            OrientedSupportKey::Plane {
                normal: a,
                offset_mm: oa,
                outward_positive: sa,
            },
            OrientedSupportKey::Plane {
                normal: b,
                offset_mm: ob,
                outward_positive: sb,
            },
        ) => oa == ob && sa == sb && a.iter().zip(b).all(|(x, y)| x.abs_diff(*y) <= 1),
        _ => false,
    }
}

fn fits_one_direction_cell(a: &[[f64; 3]], b: &[[f64; 3]]) -> bool {
    !a.is_empty()
        && !b.is_empty()
        && (0..3).all(|axis| {
            let (lo, hi) = a
                .iter()
                .chain(b)
                .map(|v| v[axis])
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                    (lo.min(v), hi.max(v))
                });
            hi - lo <= DIRECTION_Q
        })
}

/// No nearest-neighbour or partial matching: each distinct carrier key
/// has one unambiguous partner, then every connected patch and chain must
/// match, including multiplicities, endpoints and orientation.
pub(super) fn compare(
    a: &AnalyticBoundarySignature,
    b: &AnalyticBoundarySignature,
    raw_a: impl FnOnce() -> Option<Normals>,
    raw_b: impl FnOnce() -> Option<Normals>,
) -> Option<bool> {
    let keys_a = a.supports.iter().collect::<BTreeSet<_>>();
    let keys_b = b.supports.iter().collect::<BTreeSet<_>>();
    if keys_a.len() != keys_b.len() {
        return Some(false);
    }
    let missing = keys_a.difference(&keys_b).copied().collect::<Vec<_>>();
    let extra = keys_b.difference(&keys_a).copied().collect::<Vec<_>>();
    if missing.is_empty() {
        return Some(false);
    }
    let raw_a = raw_a()?;
    let raw_b = raw_b()?;
    let mut mapping = BTreeMap::new();
    let mut used = BTreeSet::new();
    for key in extra {
        let candidates = missing
            .iter()
            .copied()
            .filter(|target| {
                adjacent_plane_keys(key, target)
                    && raw_a
                        .get(*target)
                        .zip(raw_b.get(key))
                        .is_some_and(|(x, y)| fits_one_direction_cell(x, y))
            })
            .collect::<Vec<_>>();
        let [target] = candidates.as_slice() else {
            return Some(false);
        };
        if !used.insert(*target) {
            return Some(false);
        }
        mapping.insert(key.clone(), (*target).clone());
    }
    let remap = |key: &mut OrientedSupportKey| {
        if let Some(target) = mapping.get(key) {
            *key = target.clone();
        }
    };
    let mut remapped = b.clone();
    for key in &mut remapped.supports {
        remap(key);
    }
    for chain in &mut remapped.chains {
        let supports = match chain {
            BoundaryChainKey::Open { supports, .. }
            | BoundaryChainKey::Circular { supports, .. }
            | BoundaryChainKey::Spline { supports, .. } => supports,
        };
        for key in supports.iter_mut() {
            remap(key);
        }
        supports.sort();
    }
    remapped.supports.sort();
    remapped.chains.sort();
    Some(a == &remapped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacent_direction_cells_do_not_expand_the_quantum() {
        let make_key = |x| OrientedSupportKey::Plane {
            normal: [x, 1, 0],
            offset_mm: 42,
            outward_positive: true,
        };
        let a = make_key(10);
        let b = make_key(11);
        let signature = |key: &OrientedSupportKey| AnalyticBoundarySignature {
            supports: vec![key.clone()],
            chains: vec![BoundaryChainKey::Open {
                supports: [key.clone(), key.clone()],
                endpoints_mm: [[0, 0, 0], [10, 0, 0]],
            }],
        };
        let sa = signature(&a);
        let sb = signature(&b);
        let raw = |key: &OrientedSupportKey, x: f64| {
            Some(Normals::from([(
                key.clone(),
                vec![[x * DIRECTION_Q, DIRECTION_Q, 0.0]],
            )]))
        };
        assert_eq!(
            compare(&sa, &sb, || raw(&a, 10.49), || raw(&b, 10.51)),
            Some(true)
        );
        assert_eq!(
            compare(&sa, &sb, || raw(&a, 9.6), || raw(&b, 11.4)),
            Some(false)
        );
        let mut changed = sb;
        changed.chains.clear();
        assert_eq!(
            compare(&sa, &changed, || raw(&a, 10.49), || raw(&b, 10.51)),
            Some(false)
        );
    }
}
