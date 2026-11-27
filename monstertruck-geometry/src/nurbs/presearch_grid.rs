//! Shared exact-grid seed search for polynomial and rational tensor-product surfaces.
//!
//! Basis windows are evaluated once per axis, and immutable spatial samples
//! are cached by complete geometric value, mapping kind and parameter range.
//! This changes neither the sampled grid nor the subsequent Newton solver.
#![expect(
    clippy::similar_names,
    reason = "paired U/V ranges and basis windows follow tensor-product notation"
)]
use super::{
    BasisWindow, Bounded, BsplineSurface, ControlPoint, MetricSpace, ParametricSurface,
    SearchParameterHint2D, SurfaceParameterRange, Zero, algo, presearch_trace,
    seed_index::SeedSpatialIndex,
};
use std::{
    any::Any,
    cell::{Cell, RefCell},
};

/// Immutable surface samples that can be shared across inverse searches.
///
/// The complete control net and knot vectors are retained for value equality,
/// not merely a geometry address: a surface edited in place must rebuild the
/// grid. The lookup bounds each thread to eight entries.
pub(super) struct CachedPresearchGrid<C, P, const RATIONAL: bool> {
    pub(super) surface: BsplineSurface<C>,
    pub(super) range: SurfaceParameterRange,
    division: usize,
    samples: Vec<(f64, f64, P)>,
    seen: Cell<usize>,
    pub(super) index: RefCell<Option<SeedSpatialIndex>>,
}

pub(super) const MAX_PRESEARCH_GRIDS_PER_THREAD: usize = 8;
/// Defer building the index until its cost is likely to be recovered by reuse.
pub(super) const SPATIAL_INDEX_MIN_QUERIES: usize = 128;

thread_local! {
    pub(super) static PRESEARCH_GRIDS: RefCell<Vec<Box<dyn Any>>> = const { RefCell::new(Vec::new()) };
}
#[cfg(test)]
thread_local! {
    pub(super) static PRESEARCH_BUILDS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

const fn identical_parameter_ranges(
    lhs: SurfaceParameterRange,
    rhs: SurfaceParameterRange,
) -> bool {
    lhs.0.0.to_bits() == rhs.0.0.to_bits()
        && lhs.0.1.to_bits() == rhs.0.1.to_bits()
        && lhs.1.0.to_bits() == rhs.1.0.to_bits()
        && lhs.1.1.to_bits() == rhs.1.1.to_bits()
}

/// Resolve a search hint without scanning when a parameter is already known.
pub(super) fn resolve_hint(
    hint: SearchParameterHint2D,
    full_range: impl FnOnce() -> SurfaceParameterRange,
    search: impl FnOnce(SurfaceParameterRange) -> (f64, f64),
) -> (f64, f64) {
    match hint {
        SearchParameterHint2D::Parameter(u, v) => (u, v),
        SearchParameterHint2D::Range(u, v) => search((u, v)),
        SearchParameterHint2D::None => search(full_range()),
    }
}

/// Return the same discrete seed as the uniform generic scan. The const
/// mapping kind distinguishes polynomial coordinates from homogeneous ones
/// even when their underlying control-point types happen to coincide.
pub(super) fn presearch<S, C, const RATIONAL: bool>(
    surface: &S,
    bsp: &BsplineSurface<C>,
    point: S::Point,
    (urange, vrange): SurfaceParameterRange,
    division: usize,
    convert: impl Fn(C) -> S::Point,
) -> (f64, f64)
where
    S: ParametricSurface,
    C: ControlPoint<f64> + PartialEq + 'static,
    S::Point: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + 'static,
{
    let started = presearch_trace::begin();
    let (degree0, degree1) = bsp.degrees();
    let (knot_u, knot_v) = &bsp.knot_vecs;
    // The separable kernel reproduces only the non-degenerate branch of
    // `derivative_mn`; defer any zero-range surface to the generic scan so
    // the degenerate corner cases stay byte-for-byte unchanged.
    if division == 0 || knot_u.range_length().is_zero() || knot_v.range_length().is_zero() {
        return algo::surface::presearch(surface, point, (urange, vrange), division);
    }
    // Spec 014 W3: the separable fast path scans the SAME grid as
    // `algo::surface::presearch`, so it charges the same nodes -- otherwise
    // the unit would read lower purely because a surface took the fast path.
    // (The degenerate early return above delegates and is charged there.)
    algo::surface::charge_presearch_nodes(algo::surface::presearch_nodes(division));
    // Preserve the original u-major/v-minor traversal and strict-<
    // first-minimum rule, including ties and NaNs.
    let nearest = |samples: &[(f64, f64, S::Point)]| {
        let mut res = (0.0, 0.0);
        let mut min = f64::INFINITY;
        for &(u, v, position) in samples {
            let dist = position.distance2(point);
            if dist < min {
                min = dist;
                res = (u, v);
            }
        }
        res
    };
    if let Some((result, indexed, built)) = PRESEARCH_GRIDS.with(|grids| {
        grids.borrow().iter().rev().find_map(|entry| {
            let entry = entry.downcast_ref::<CachedPresearchGrid<C, S::Point, RATIONAL>>()?;
            (entry.division == division
                && identical_parameter_ranges(entry.range, (urange, vrange))
                && entry.surface == *bsp)
                .then(|| {
                    let count = entry.seen.get().saturating_add(1);
                    entry.seen.set(count);
                    let built = count == SPATIAL_INDEX_MIN_QUERIES;
                    if built {
                        *entry.index.borrow_mut() = SeedSpatialIndex::build(&entry.samples);
                    }
                    let index = entry
                        .index
                        .borrow()
                        .as_ref()
                        .and_then(|tree| tree.nearest(&entry.samples, point));
                    let indexed = index.is_some();
                    let result = index.map_or_else(
                        || nearest(&entry.samples),
                        |i| (entry.samples[i].0, entry.samples[i].1),
                    );
                    (result, indexed, built)
                })
        })
    }) {
        presearch_trace::seed(true, indexed, built, started);
        return result;
    }

    let ((u0, u1), (v0, v1)) = (urange, vrange);
    let vbasis: Vec<BasisWindow> = (0..=division)
        .map(|j| {
            let q = crate::usize_to_f64(j) / crate::usize_to_f64(division);
            let v = f64::mul_add(v1, q, v0 * (1.0 - q));
            knot_v.bspline_basis_functions_unchecked(degree1, 0, v)
        })
        .collect();
    let mut samples = Vec::with_capacity(vbasis.len().saturating_mul(vbasis.len()));
    for i in 0..=division {
        let p = crate::usize_to_f64(i) / crate::usize_to_f64(division);
        let u = f64::mul_add(u1, p, u0 * (1.0 - p));
        let ubasis = knot_u.bspline_basis_functions_unchecked(degree0, 0, u);
        for (j, vb) in vbasis.iter().enumerate() {
            let q = crate::usize_to_f64(j) / crate::usize_to_f64(division);
            let v = f64::mul_add(v1, q, v0 * (1.0 - q));
            let position = convert(C::from_vec(bsp.combine_basis(&ubasis, vb)));
            samples.push((u, v, position));
        }
    }
    let result = nearest(&samples);
    #[cfg(test)]
    PRESEARCH_BUILDS.with(|count| count.set(count.get() + 1));
    PRESEARCH_GRIDS.with(|grids| {
        let mut grids = grids.borrow_mut();
        if grids.len() == MAX_PRESEARCH_GRIDS_PER_THREAD {
            grids.remove(0);
        }
        grids.push(Box::new(CachedPresearchGrid::<C, S::Point, RATIONAL> {
            surface: bsp.clone(),
            range: (urange, vrange),
            division,
            samples,
            seen: Cell::new(1),
            index: RefCell::new(None),
        }));
    });
    presearch_trace::seed(false, false, false, started);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::{BoundedSurface, KnotVector, Point3};

    fn polynomial_surface() -> BsplineSurface<Point3> {
        BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(3), KnotVector::bezier_knot(1)),
            vec![
                vec![Point3::new(2., 0., 0.), Point3::new(2., 0., 1.)],
                vec![Point3::new(2.5, 1., 0.), Point3::new(2.5, 1., 1.)],
                vec![Point3::new(3.5, 0.2, 0.), Point3::new(3.5, 0.2, 1.)],
                vec![Point3::new(4., 1., 0.), Point3::new(4., 1., 1.)],
            ],
        )
    }

    fn check(
        surface: &BsplineSurface<Point3>,
        point: Point3,
        range: SurfaceParameterRange,
        division: usize,
    ) {
        let actual = presearch::<_, _, false>(
            surface,
            surface,
            point,
            range,
            division,
            std::convert::identity,
        );
        let expected = algo::surface::presearch(surface, point, range, division);
        assert_eq!(
            [actual.0.to_bits(), actual.1.to_bits()],
            [expected.0.to_bits(), expected.1.to_bits()]
        );
    }

    #[test]
    fn polynomial_presearch_reuses_exact_grid_and_index() {
        let surface = polynomial_surface();
        let full = surface.range_tuple();
        let point = Point3::new(2.6, 0.3, 0.7);
        let builds = || PRESEARCH_BUILDS.with(Cell::get);
        let initial = builds();
        for i in 0..SPATIAL_INDEX_MIN_QUERIES + 2 {
            let q = Point3::new(2.1 + crate::usize_to_f64(i % 61) * 0.03, 0.25, 0.5);
            check(&surface, q, full, 50);
        }
        assert_eq!(builds(), initial + 1);
        PRESEARCH_GRIDS.with(|grids| {
            let grids = grids.borrow();
            let entry = grids
                .iter()
                .filter_map(|entry| {
                    entry.downcast_ref::<CachedPresearchGrid<Point3, Point3, false>>()
                })
                .find(|entry| entry.surface == surface && entry.range == full)
                .expect("polynomial grid");
            assert!(entry.index.borrow().is_some());
        });
        // Nonfinite queries retain the original scan's default result.
        check(&surface, Point3::new(f64::NAN, 0., 0.), full, 50);
        // Range, resolution and geometry are independent cache-key dimensions.
        check(&surface, point, ((0.2, 0.7), (0.1, 0.6)), 50);
        check(&surface, point, full, 12);
        let mut changed = surface.clone();
        changed.control_points[1][0].y += 0.07;
        check(&changed, point, full, 50);
        assert_eq!(builds(), initial + 4);
        check(&changed, point, full, 0);
        assert_eq!(builds(), initial + 4);
    }
}
