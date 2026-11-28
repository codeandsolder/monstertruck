//! Canonical cells are a representation, not geometric discontinuities.
//! Only independently certified signatures may be compared here. A cell
//! crossing requires unrounded witnesses within the physical recovery
//! budget, unique correspondence and complete signature equality.
use super::{
    AnalyticBoundarySignature, BoundaryChainKey, DIRECTION_Q, OrientedSupportKey, POINT_Q_MM,
    PointKey, SCALAR_Q_MM, canonical_boundary_direction, point_key, quantize,
};
use num_traits::ToPrimitive;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum CellKind {
    Point,
    Direction,
    Offset,
    Radius,
    Slope,
}

impl CellKind {
    const fn comparison_limit(self, tolerance_mm: f64) -> f64 {
        match self {
            Self::Point | Self::Offset | Self::Radius => tolerance_mm,
            Self::Direction | Self::Slope => DIRECTION_Q,
        }
    }

    const fn quantum(self) -> f64 {
        match self {
            Self::Point => POINT_Q_MM,
            Self::Direction | Self::Slope => DIRECTION_Q,
            Self::Offset | Self::Radius => SCALAR_Q_MM,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct CellKey {
    kind: CellKind,
    components: [i64; 3],
}

impl CellKey {
    fn possibly_matches(self, other: Self, tolerance_mm: f64) -> bool {
        self.kind == other.kind
            && self.components.iter().zip(other.components).all(|(a, b)| {
                a.abs_diff(b).to_f64().is_some_and(|delta| {
                    delta <= self.kind.comparison_limit(tolerance_mm) / self.kind.quantum() + 1.0
                })
            })
    }
}

type RawBounds = ([f64; 3], [f64; 3]);
type BoundsByCell = BTreeMap<CellKey, RawBounds>;

/// Track the actual values used to construct cells, not their decoded centres.
/// Ordinary equal-signature proofs do not allocate this map; it is collected
/// only when a complete geometric correspondence could be relevant.
pub(super) struct BoundaryProofContext {
    bounds: Option<BoundsByCell>,
    tolerance_mm: f64,
}

impl BoundaryProofContext {
    pub(super) const fn new(tolerance_mm: f64) -> Self {
        Self { bounds: None, tolerance_mm }
    }

    pub(super) const fn tracked(tolerance_mm: f64) -> Self {
        Self { bounds: Some(BTreeMap::new()), tolerance_mm }
    }

    pub(super) const fn tolerance_mm(&self) -> f64 { self.tolerance_mm }

    pub(super) const fn enabled(&self) -> bool { self.bounds.is_some() }

    fn record(&mut self, kind: CellKind, components: [i64; 3], raw: [f64; 3]) {
        if let Some(bounds) = &mut self.bounds {
            let (lo, hi) = bounds
                .entry(CellKey { kind, components })
                .or_insert((raw, raw));
            for i in 0..3 {
                lo[i] = lo[i].min(raw[i]);
                hi[i] = hi[i].max(raw[i]);
            }
        }
    }

    pub(super) fn record_point(&mut self, key: PointKey, point: [f64; 3]) {
        self.record(CellKind::Point, key, point);
    }

    pub(super) fn point(&mut self, point: [f64; 3]) -> Option<PointKey> {
        let key = point_key(point)?;
        self.record_point(key, point);
        Some(key)
    }

    pub(super) fn record_direction(&mut self, key: [i64; 3], direction: [f64; 3]) {
        self.record(CellKind::Direction, key, direction);
    }

    pub(super) fn direction(&mut self, raw: [f64; 3]) -> Option<([f64; 3], [i64; 3])> {
        let (direction, key) = canonical_boundary_direction(raw)?;
        self.record_direction(key, direction);
        Some((direction, key))
    }

    pub(super) fn scalar(&mut self, kind: CellKind, value: f64) -> Option<i64> {
        let key = quantize(value, kind.quantum())?;
        self.record(kind, [key, 0, 0], [value, 0., 0.]);
        Some(key)
    }

    fn fits(&self, a: CellKey, other: &Self, b: CellKey, tolerance_mm: f64) -> bool {
        let Some((alo, ahi)) = self.bounds.as_ref().and_then(|cells| cells.get(&a)) else {
            return false;
        };
        let Some((blo, bhi)) = other.bounds.as_ref().and_then(|cells| cells.get(&b)) else {
            return false;
        };
        let span: [f64; 3] = std::array::from_fn(|i| ahi[i].max(bhi[i]) - alo[i].min(blo[i]));
        match a.kind {
            CellKind::Point => crate::math3::norm(span) <= a.kind.comparison_limit(tolerance_mm),
            _ => span
                .iter()
                .all(|&d| d <= a.kind.comparison_limit(tolerance_mm)),
        }
    }
}

#[cfg(test)]
impl Default for BoundaryProofContext {
    fn default() -> Self {
        Self::new(super::TEST_BOUNDARY_TOL_MM)
    }
}

fn map_vector(key: &mut [i64; 3], kind: CellKind, map: &mut dyn FnMut(CellKey) -> CellKey) {
    *key = map(CellKey {
        kind,
        components: *key,
    })
    .components;
}
fn map_scalar(key: &mut i64, kind: CellKind, map: &mut dyn FnMut(CellKey) -> CellKey) {
    *key = map(CellKey {
        kind,
        components: [*key, 0, 0],
    })
    .components[0];
}

fn map_support(key: &mut OrientedSupportKey, map: &mut dyn FnMut(CellKey) -> CellKey) {
    use CellKind::{Direction, Offset, Point, Radius, Slope};
    match key {
        OrientedSupportKey::Plane {
            normal, offset_mm, ..
        } => {
            map_vector(normal, Direction, map);
            map_scalar(offset_mm, Offset, map);
        }
        OrientedSupportKey::Cylinder {
            axis,
            closest_mm,
            radius_mm,
            ..
        } => {
            map_vector(axis, Direction, map);
            map_vector(closest_mm, Point, map);
            map_scalar(radius_mm, Radius, map);
        }
        OrientedSupportKey::Cone {
            apex_mm,
            axis,
            slope,
            ..
        } => {
            map_vector(apex_mm, Point, map);
            map_vector(axis, Direction, map);
            map_scalar(slope, Slope, map);
        }
        OrientedSupportKey::Sphere {
            center_mm,
            radius_mm,
            ..
        } => {
            map_vector(center_mm, Point, map);
            map_scalar(radius_mm, Radius, map);
        }
        OrientedSupportKey::Torus {
            center_mm,
            axis,
            major_radius_mm,
            minor_radius_mm,
            ..
        } => {
            map_vector(center_mm, Point, map);
            map_vector(axis, Direction, map);
            map_scalar(major_radius_mm, Radius, map);
            map_scalar(minor_radius_mm, Radius, map);
        }
        OrientedSupportKey::SplineStrip { strip, .. } => {
            strip.map_points(|key| map_vector(key, Point, map));
        }
    }
}

const fn chain_supports(chain: &mut BoundaryChainKey) -> &mut [OrientedSupportKey; 2] {
    match chain {
        BoundaryChainKey::Open { supports, .. }
        | BoundaryChainKey::Circular { supports, .. }
        | BoundaryChainKey::CylinderSection { supports, .. }
        | BoundaryChainKey::CylinderIntersection { supports, .. }
        | BoundaryChainKey::Spline { supports, .. } => supports,
    }
}

fn map_chain_geometry(chain: &mut BoundaryChainKey, map: &mut dyn FnMut(CellKey) -> CellKey) {
    match chain {
        BoundaryChainKey::Open { endpoints_mm, .. } => {
            for key in endpoints_mm {
                map_vector(key, CellKind::Point, map);
            }
        }
        BoundaryChainKey::Circular {
            locus,
            endpoints_mm,
            ..
        } => {
            map_vector(&mut locus.center_mm, CellKind::Point, map);
            map_vector(&mut locus.normal, CellKind::Direction, map);
            map_scalar(&mut locus.radius_mm, CellKind::Radius, map);
            for key in endpoints_mm.iter_mut().flatten() {
                map_vector(key, CellKind::Point, map);
            }
        }
        BoundaryChainKey::CylinderSection { endpoints_mm, .. }
        | BoundaryChainKey::CylinderIntersection { endpoints_mm, .. } => {
            for key in endpoints_mm.iter_mut().flatten() {
                map_vector(key, CellKind::Point, map);
            }
        }
        BoundaryChainKey::Spline { curve, .. } => {
            curve.map_points(|key| map_vector(key, CellKind::Point, map));
        }
    }
}

fn chain_endpoints(chain: &BoundaryChainKey) -> Result<Option<[PointKey; 2]>, ()> {
    Ok(match chain {
        BoundaryChainKey::Open { endpoints_mm, .. } => Some(*endpoints_mm),
        BoundaryChainKey::Circular { endpoints_mm, .. }
        | BoundaryChainKey::CylinderSection { endpoints_mm, .. }
        | BoundaryChainKey::CylinderIntersection { endpoints_mm, .. } => *endpoints_mm,
        BoundaryChainKey::Spline { curve, .. } => Some(curve.ends().ok_or(())?),
    })
}

fn endpoints(signature: &AnalyticBoundarySignature) -> Option<BTreeSet<PointKey>> {
    let mut points = BTreeSet::new();
    for chain in &signature.chains {
        points.extend(chain_endpoints(chain).ok()?.into_iter().flatten());
    }
    Some(points)
}

/// Check corresponding typed fields in structural order, then require exact
/// equality of the complete record. Rational weights, orientation, nappe,
/// span structure and all other non-cell fields are never rewritten.
fn structural_match<T: Clone + Eq>(
    a: &T,
    b: &T,
    walk: impl Fn(&mut T, &mut dyn FnMut(CellKey) -> CellKey),
    fits: &impl Fn(CellKey, CellKey) -> bool,
) -> bool {
    let mut fields = Vec::new();
    walk(&mut a.clone(), &mut |key| {
        fields.push(key);
        key
    });
    let mut fields = fields.into_iter();
    let mut remapped = b.clone();
    let mut matched = true;
    walk(&mut remapped, &mut |key| {
        if let Some(target) = fields.next()
            && fits(target, key)
        {
            target
        } else {
            matched = false;
            key
        }
    });
    matched && fields.next().is_none() && a == &remapped
}

fn chain_matches(
    a: &BoundaryChainKey,
    b: &BoundaryChainKey,
    point_map: &BTreeMap<PointKey, PointKey>,
    fits: &impl Fn(CellKey, CellKey) -> bool,
) -> bool {
    let Ok(ends_a) = chain_endpoints(a) else {
        return false;
    };
    let Ok(ends_b) = chain_endpoints(b) else {
        return false;
    };
    let mapped_ends =
        ends_b.map(|ends| ends.map(|point| point_map.get(&point).copied().unwrap_or(point)));
    if ends_a == mapped_ends && structural_match(a, b, map_chain_geometry, fits) {
        return true;
    }
    // Linear chains have unordered endpoints; circular and spline chains
    // retain the directed canonical interval and must not be reversed.
    let mut reversed = b.clone();
    if let BoundaryChainKey::Open { endpoints_mm, .. } = &mut reversed {
        endpoints_mm.reverse();
        ends_a == mapped_ends.map(|[first, second]| [second, first])
            && structural_match(a, &reversed, map_chain_geometry, fits)
    } else {
        false
    }
}

fn counts<T: Ord>(items: &[T]) -> BTreeMap<&T, usize> {
    let mut counts = BTreeMap::new();
    for key in items {
        *counts.entry(key).or_insert(0) += 1;
    }
    counts
}

/// Fix exact identities, then require a unique bijection of distinct complete
/// records with unchanged multiplicities. Direction/radius/offset values
/// themselves need not be unique: two different carriers can share a value.
fn correspondence<T: Ord + Clone>(
    a: &[T],
    b: &[T],
    fits: impl Fn(&T, &T) -> bool,
) -> Option<BTreeMap<T, T>> {
    let counts_a = counts(a);
    let counts_b = counts(b);
    let mut missing = BTreeMap::new();
    for (&key, &count) in &counts_a {
        match counts_b.get(key) {
            Some(&other) if count == other => {}
            Some(_) => return None,
            None => {
                missing.insert(key, count);
            }
        }
    }
    let mut mapping = BTreeMap::new();
    let mut used = BTreeSet::new();
    for (&key, &count) in &counts_b {
        if counts_a.contains_key(key) {
            continue;
        }
        let mut candidates = missing
            .iter()
            .filter(|(target, n)| **n == count && fits(target, key));
        let (&target, _) = candidates.next()?;
        if candidates.next().is_some() || !used.insert(target) {
            return None;
        }
        mapping.insert(key.clone(), target.clone());
    }
    (used.len() == missing.len()).then_some(mapping)
}

const fn point_cell(components: PointKey) -> CellKey {
    CellKey {
        kind: CellKind::Point,
        components,
    }
}

pub(super) fn compare(
    a: &AnalyticBoundarySignature,
    b: &AnalyticBoundarySignature,
    tolerance_mm: f64,
    raw_a: impl FnOnce() -> Option<BoundaryProofContext>,
    raw_b: impl FnOnce() -> Option<BoundaryProofContext>,
) -> Option<bool> {
    if a.supports.len() != b.supports.len() || a.chains.len() != b.chains.len() {
        return Some(false);
    }
    let points_a = endpoints(a)?;
    let points_b = endpoints(b)?;
    // Cheap negative evidence before re-running geometric proofs. A coarse
    // ambiguity is not a rejection: raw intervals may disambiguate it.
    if points_a.len() != points_b.len() {
        return Some(false);
    }
    let missing = points_a.difference(&points_b).copied().collect::<Vec<_>>();
    if points_b.difference(&points_a).any(|point| {
        !missing
            .iter()
            .any(|target| point_cell(*target).possibly_matches(point_cell(*point), tolerance_mm))
    }) {
        return Some(false);
    }
    let points_a = points_a.into_iter().collect::<Vec<_>>();
    let points_b = points_b.into_iter().collect::<Vec<_>>();
    let raw_a = raw_a()?;
    let raw_b = raw_b()?;
    let raw_fits = |target: CellKey, key: CellKey| {
        target == key
            || (target.possibly_matches(key, tolerance_mm)
                && raw_a.fits(target, &raw_b, key, tolerance_mm))
    };
    let point_map = correspondence(&points_a, &points_b, |target, key| {
        raw_fits(point_cell(*target), point_cell(*key))
    })?;
    let support_map = correspondence(&a.supports, &b.supports, |target, key| {
        structural_match(target, key, map_support, &raw_fits)
    })?;
    let mut remapped = b.chains.clone();
    for chain in &mut remapped {
        let supports = chain_supports(chain);
        for key in supports.iter_mut() {
            if let Some(target) = support_map.get(key) {
                *key = target.clone();
            }
        }
        let reordered = supports[0] > supports[1];
        supports.sort();
        if reordered
            && let BoundaryChainKey::CylinderIntersection {
                primary_is_second, ..
            } = chain
        {
            *primary_is_second = !*primary_is_second;
        }
    }
    Some(
        correspondence(&a.chains, &remapped, |target, key| {
            chain_matches(target, key, &point_map, &raw_fits)
        })
        .is_some(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane_signature(
        normal: [i64; 3],
        offset: i64,
        endpoint: PointKey,
    ) -> AnalyticBoundarySignature {
        let key = OrientedSupportKey::Plane {
            normal,
            offset_mm: offset,
            outward_positive: true,
        };
        AnalyticBoundarySignature {
            supports: vec![key.clone()],
            chains: vec![BoundaryChainKey::Open {
                supports: [key.clone(), key],
                endpoints_mm: [[0, 0, 0], endpoint],
            }],
        }
    }

    #[test]
    fn raw_direction_offset_and_point_cells_match_together() {
        let a = plane_signature([10, 1, 0], 20, [30, 0, 0]);
        let b = plane_signature([11, 1, 0], 21, [31, 0, 0]);
        let raw = |delta: f64| {
            let mut cells = BoundaryProofContext::tracked(POINT_Q_MM);
            cells.record(
                CellKind::Direction,
                [quantize((10. + delta) * DIRECTION_Q, DIRECTION_Q)?, 1, 0],
                [(10. + delta) * DIRECTION_Q, DIRECTION_Q, 0.],
            );
            cells.scalar(CellKind::Offset, (20. + delta) * SCALAR_Q_MM)?;
            cells.point([(30. + delta) * POINT_Q_MM, 0., 0.])?;
            Some(cells)
        };
        assert_eq!(
            compare(&a, &b, POINT_Q_MM, || raw(0.49), || raw(0.51)),
            Some(true)
        );
        assert_ne!(
            compare(&a, &b, POINT_Q_MM, || raw(-0.4), || raw(1.4)),
            Some(true)
        );
        let mut changed = b.clone();
        changed.chains.clear();
        assert_eq!(
            compare(&a, &changed, POINT_Q_MM, || raw(0.49), || raw(0.51)),
            Some(false)
        );
        assert_ne!(
            compare(
                &a,
                &b,
                POINT_Q_MM,
                || Some(BoundaryProofContext::tracked(POINT_Q_MM)),
                || raw(0.51)
            ),
            Some(true)
        );
    }

    #[test]
    fn ambiguous_and_colliding_cell_matches_are_rejected() {
        let a = [0_i64, 2];
        let b = [1_i64, 3];
        assert!(correspondence(&a, &b, |target, key| target.abs_diff(*key) <= 1).is_none());
        assert!(correspondence(&a, &b, |target, _| *target == 2).is_none());
        assert!(correspondence(&a, &b, |target, key| target + 1 == *key).is_some());
    }

    #[test]
    fn distinct_carriers_may_share_direction_cells() {
        let mut a = plane_signature([10, 1, 0], 20, [30, 0, 0]);
        let mut b = plane_signature([11, 1, 0], 20, [30, 0, 0]);
        let common = plane_signature([11, 1, 0], 50, [60, 0, 0]);
        for signature in [&mut a, &mut b] {
            signature.supports.extend(common.supports.clone());
            signature.chains.extend(common.chains.clone());
        }
        let raw = |key, x: f64| {
            let mut cells = BoundaryProofContext::tracked(POINT_Q_MM);
            cells.record(
                CellKind::Direction,
                [key, 1, 0],
                [x * DIRECTION_Q, DIRECTION_Q, 0.],
            );
            Some(cells)
        };
        assert_eq!(
            compare(&a, &b, POINT_Q_MM, || raw(10, 10.49), || raw(11, 10.51)),
            Some(true)
        );
    }
}
