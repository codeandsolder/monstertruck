//! Face boundaries in surface parameter space.
//!
//! The work splits four ways, one file each: [`projection`] turns a face's 3D
//! wire into a uv polyline, [`projection_debug`] is the measurement-only lens
//! over that step, [`periodic`] normalises loops on periodic and singular
//! domains, and [`loops`] assembles the pieces into the boundary the
//! triangulator consumes. What stays here is the vocabulary all four share.

use super::{
    BoundedCurve, Cdt, EuclideanSpace, ExactParameterBoundary2D, ExactTrimBoundary2D,
    FixedVertexHandle, HashGen, HashMap, InnerSpace, IsoparametricCurveOptions, Itertools,
    KERNEL_EPSILON, MetricSpace, Origin, ParameterBoundary2D, ParametricCurve3D,
    ParametricSurface3D, Point2, Point3, PolylineCurve, PreMeshableSurface, SP, SPoint2, Tolerance,
    Triangulation, Vector2, algo, env, mesh, resample_boundary, simplify_parameter_boundary,
};

mod loops;
mod periodic;
mod projection;
mod projection_debug;

// `periodic` is used only within `boundary`, so its glob stays module-private.
use periodic::{
    align_full_period_open_pair, align_periodic_open_pair, close_open_periodic_curve_to_singular,
    compact_periodic_closed_loop, normalize_range, open_periodic_closed_loop,
    periodic_axis_full_span, periodic_min_difference, periodic_open_axis_full_span,
};

pub(super) use loops::*;
pub(super) use projection::*;
pub(super) use projection_debug::*;

#[derive(Clone, Copy, Debug, derive_more::Deref, derive_more::DerefMut)]
pub(super) struct SurfacePoint {
    pub(super) point: Point3,
    #[deref]
    #[deref_mut]
    pub(super) uv: Point2,
}

impl From<(Point2, Point3)> for SurfacePoint {
    fn from((uv, point): (Point2, Point3)) -> Self { Self { point, uv } }
}

fn connect_edges<P>(vecs: impl IntoIterator<Item = Vec<P>>) -> Vec<P> {
    let closure = |vec: Vec<P>| {
        let len = vec.len();
        vec.into_iter().take(len - 1)
    };
    vecs.into_iter().flat_map(closure).collect()
}

pub(super) type UvKey = (u64, u64);

pub(super) const fn uv_key(uv: Point2) -> UvKey { (uv.x.to_bits(), uv.y.to_bits()) }

pub(super) fn surface_point_with_cache(
    surface: &impl PreMeshableSurface,
    uv: Point2,
    point_cache: &mut HashMap<UvKey, Point3>,
) -> SurfacePoint {
    let point = *point_cache
        .entry(uv_key(uv))
        .or_insert_with(|| surface.evaluate(uv.x, uv.y));
    (uv, point).into()
}

fn polyline_on_surface(
    surface: impl PreMeshableSurface,
    p: SurfacePoint,
    q: SurfacePoint,
    tolerance: f64,
    point_cache: &mut HashMap<UvKey, Point3>,
) -> Vec<SurfacePoint> {
    use monstertruck_geometry::prelude::*;
    let line = Line(p.uv, q.uv);
    let pcurve = ParameterCurve::new(line, &surface);
    let (vec, _) = pcurve.parameter_division(pcurve.range_tuple(), tolerance);
    vec.into_iter()
        .map(|t| {
            let uv = line.evaluate(t);
            surface_point_with_cache(&surface, uv, point_cache)
        })
        .collect()
}

#[cfg(test)]
mod tests;
