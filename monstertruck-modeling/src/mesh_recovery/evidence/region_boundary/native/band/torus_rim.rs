//! Variable rims on regular toroidal bands. Lifted source chords are
//! tensor products of two short, directed rational circle arcs; their
//! provenance supplies the continuous chart path, not just endpoint samples.
//! A linear shear of the unwrapped coordinates can certify a simple cycle
//! even when one edge advances slightly backwards in the periodic direction.

use super::{BandContext, BandPlan, Ring, copy_rings, parameter_circle};
use crate::mesh_recovery::{OrientedTrim, SharedTrimGeometry, SurfaceCarrier};
use crate::{BoundedCurve, Edge, NurbsCurve, ParametricCurve, Point3};
use monstertruck_geometry::analytic_curve_incidence::{
    CircularArcFrame, certify_circular_arc_coverage,
};
use monstertruck_geometry::prelude::TryIntoHomogeneousBsplineCurve;
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy, Debug)]
enum PathKind {
    Lifted,
    Circle,
}

#[derive(Clone, Copy, Debug)]
struct RimPath {
    advance: f64,
    height: [f64; 2],
    kind: PathKind,
}

fn signed_delta(a: f64, b: f64) -> f64 { (b - a + PI).rem_euclid(TAU) - PI }

// For the symmetric positive-weight quadratic unit-circle arc of signed
// sweep d, angular speed is 2*sin(d/2)/w(t), with w(t) ranging between
// (1+cos(d/2))/2 and 1. Both factors in a lifted torus chord use the same t.
fn speed_bounds(delta: f64) -> [f64; 2] {
    let a = 2.0 * (delta * 0.5).sin();
    let b = 4.0 * (delta * 0.25).tan();
    [a.min(b), a.max(b)]
}

fn shear_feasible(paths: &[RimPath], sign: f64, resolution: f64) -> bool {
    let mut lower = 0.0_f64;
    let mut upper = f64::INFINITY;
    // For circle spans only one coordinate varies, so its directed
    // monotonicity is independent of its rational parameter speed.
    for path in paths {
        let dh = path.height[1] - path.height[0];
        let (a, b) = match path.kind {
            PathKind::Lifted => (
                speed_bounds(path.advance)[0] - resolution,
                speed_bounds(sign * dh)[0] - resolution,
            ),
            PathKind::Circle if path.advance > resolution => continue,
            PathKind::Circle if path.advance == 0.0 && sign * dh > resolution => continue,
            PathKind::Circle => return false,
        };
        // Require d(u + sign*k*h)/dt > resolution*(1+k) on every
        // lifted span. Intersect the open linear inequalities in k.
        if b > 0.0 {
            lower = lower.max(-a / b);
        } else if b < 0.0 {
            upper = upper.min(-a / b);
        } else if a <= 0.0 {
            return false;
        }
    }
    lower.is_finite() && lower < upper && upper > 0.0
}

fn simple_cycle(paths: &[RimPath], resolution: f64) -> bool {
    let sweep: f64 = paths.iter().map(|p| p.advance).sum();
    if (sweep - TAU).abs() > resolution {
        return false;
    }
    // A strictly progressing unwrapped angular coordinate works directly.
    // Constant-coordinate circular steps are handled by a nonzero shear.
    paths.iter().all(|p| match p.kind {
        PathKind::Lifted => speed_bounds(p.advance)[0] > resolution,
        PathKind::Circle => p.advance > resolution,
    }) || shear_feasible(paths, 1.0, resolution)
        || shear_feasible(paths, -1.0, resolution)
}

// A full meridional winding requires a ring torus. An azimuthal
// annulus may also occupy a regular positive-radius part of a horn or
// spindle torus, provided the entire meridian interval excludes poles.
fn regular_radial_band(
    major: f64,
    minor: f64,
    varying: usize,
    lower: f64,
    span: f64,
    tol: f64,
) -> bool {
    if major <= crate::KERNEL_EPSILON || minor <= tol {
        return false;
    }
    let minimum = if varying == 1 || (PI - lower).rem_euclid(TAU) <= span {
        -1.0
    } else {
        lower.cos().min((lower + span).cos())
    };
    minor.mul_add(minimum, major) > tol
}

impl BandContext<'_> {
    fn torus_chart(&self, point: Point3) -> Option<[f64; 2]> {
        super::chart_point(self.carrier, point, self.options.angular_resolution)
    }

    fn circle_delta(&self, edge: &Edge, uv: [f64; 2]) -> Option<[f64; 2]> {
        let curve = NurbsCurve::new(edge.oriented_curve().try_into_homogeneous_bspline_curve()?);
        let tolerance = self.options.intersection_tolerance;
        for axis in 0..2 {
            let Some(circle) = parameter_circle(self.carrier, axis, uv[1 - axis]) else {
                continue;
            };
            let Some(frame) = CircularArcFrame::new(
                circle.center,
                circle.normal,
                circle.radius,
                edge.front().point(),
                tolerance,
            ) else {
                continue;
            };
            let Some(coverage) = certify_circular_arc_coverage(&curve, frame) else {
                continue;
            };
            let angular = self
                .options
                .angular_resolution
                .max(tolerance / circle.radius);
            let positive = coverage.start_rad <= angular;
            let negative = (coverage.start_rad + coverage.sweep_rad - TAU).abs() <= angular;
            if positive == negative {
                continue;
            }
            let mut delta = [0.0; 2];
            delta[axis] = if positive {
                coverage.sweep_rad
            } else {
                -coverage.sweep_rad
            };
            return Some(delta);
        }
        None
    }

    fn torus_path(&self, edge: &Edge, kind: PathKind, ring: Ring, side: f64) -> Option<RimPath> {
        let start = self.torus_chart(edge.front().point())?;
        let end = self.torus_chart(edge.back().point())?;
        let resolution = self.options.angular_resolution;
        let delta = match kind {
            PathKind::Lifted => {
                let delta = std::array::from_fn(|i| signed_delta(start[i], end[i]));
                if delta.iter().any(|d| d.abs() >= PI - resolution) {
                    return None;
                }
                delta
            }
            PathKind::Circle => self.circle_delta(edge, start)?,
        };
        if (0..2).any(|i| signed_delta(start[i] + delta[i], end[i]).abs() > resolution) {
            return None;
        }
        let fixed = 1 - ring.axis;
        let h0 = (side * (start[fixed] - ring.fixed)).rem_euclid(TAU);
        let h1 = side.mul_add(delta[fixed], h0);
        let height = [h0, h1];
        if height
            .iter()
            .any(|h| !h.is_finite() || *h <= resolution || *h >= TAU - resolution)
        {
            return None;
        }
        let direction = if ring.positive { -1.0 } else { 1.0 };
        Some(RimPath {
            advance: direction * delta[ring.axis],
            height,
            kind,
        })
    }

    fn torus_cuts(
        &self,
        rim: &[Edge],
        paths: &[RimPath],
        ring: Ring,
    ) -> Option<[([f64; 2], Point3); 2]> {
        let resolution = self.options.angular_resolution;
        let direction = if ring.positive { -1.0 } else { 1.0 };
        let first = self.torus_chart(rim.first()?.front().point())?[ring.axis];
        let mut position = 0.0;
        let ranges = paths
            .iter()
            .map(|p| {
                let range = [position, position + p.advance];
                position += p.advance;
                range
            })
            .collect::<Vec<_>>();
        let mut candidates = Vec::new();
        for (i, edge) in rim.iter().enumerate() {
            if paths[i].advance <= 2.0 * resolution {
                continue;
            }
            let curve = edge.oriented_curve();
            let (a, b) = curve.range_tuple();
            let point = curve.evaluate(a.midpoint(b));
            let uv = self.torus_chart(point)?;
            let phase = (direction * (uv[ring.axis] - first)).rem_euclid(TAU);
            let [lo, hi] = ranges[i];
            let angle = ((lo.midpoint(hi) - phase) / TAU)
                .round()
                .mul_add(TAU, phase);
            if angle <= lo + resolution || angle >= hi - resolution {
                continue;
            }
            // The chosen iso-coordinate intersects the *entire* second
            // rim just once, including any backwards-going steps.
            if ranges.iter().enumerate().any(|(j, &[a, b])| {
                j != i && {
                    let p = ((a.midpoint(b) - angle) / TAU).round().mul_add(TAU, angle);
                    p >= a.min(b) - resolution && p <= a.max(b) + resolution
                }
            }) {
                continue;
            }
            candidates.push((uv, point));
        }
        let first = *candidates.first()?;
        let distance =
            |p: &([f64; 2], Point3)| signed_delta(first.0[ring.axis], p.0[ring.axis]).abs();
        let second = *candidates
            .iter()
            .max_by(|a, b| distance(a).total_cmp(&distance(b)))?;
        (distance(&second) > resolution).then_some([first, second])
    }

    pub(super) fn variable_torus_band(
        &self,
        cycles: [&[OrientedTrim]; 2],
        ring_index: usize,
        ring: Ring,
        trims: &mut [Result<Vec<Edge>, super::NativeTrimError>],
    ) -> Option<BandPlan> {
        let SurfaceCarrier::CircularRevolution {
            radial_center,
            profile_radius,
            ..
        } = self.carrier
        else {
            return None;
        };
        if radial_center <= crate::KERNEL_EPSILON {
            return None; // Spherical chord lifts are not torus tensor-product arcs.
        }
        let changed = copy_rings(cycles[0].iter().chain(cycles[1]), trims)?;
        let side = if (ring.positive == (ring.axis == 0)) == self.outward {
            1.0
        } else {
            -1.0
        };
        let mut rim = Vec::new();
        let mut paths = Vec::new();
        for reference in cycles[1 - ring_index] {
            let edges = changed.get(&reference.trim)?;
            let trim = self.graph.trims.get(reference.trim)?;
            let kind = if self
                .lifted_deviations
                .get(reference.trim)
                .is_some_and(Option::is_some)
            {
                if edges.len() != trim.source_edges.len() {
                    return None;
                }
                PathKind::Lifted
            } else if matches!(trim.geometry, SharedTrimGeometry::Circle(_)) {
                PathKind::Circle
            } else {
                return None;
            };
            let oriented: Vec<_> = if reference.reversed {
                edges.iter().rev().map(Edge::inverse).collect()
            } else {
                edges.clone()
            };
            for edge in oriented {
                paths.push(self.torus_path(&edge, kind, ring, side)?);
                rim.push(edge);
            }
        }
        let resolution = self.options.angular_resolution;
        if !simple_cycle(&paths, resolution)
            || paths
                .iter()
                .zip(paths.iter().cycle().skip(1))
                .any(|(a, b)| (a.height[1] - b.height[0]).abs() > resolution)
        {
            return None;
        }
        let extent = paths.iter().flat_map(|p| p.height).fold(0.0_f64, f64::max);
        let lower = if side > 0.0 {
            ring.fixed
        } else {
            ring.fixed - extent
        };
        if !regular_radial_band(
            radial_center,
            profile_radius,
            ring.axis,
            lower,
            extent,
            self.options.intersection_tolerance,
        ) || !self.covers_band(1 - ring.axis, lower, extent)
            || !self.supports.iter().all(|&support| {
                self.evidence.planar_supports[support]
                    .source_face_indices
                    .iter()
                    .all(|&polygon| {
                        self.evidence
                            .source_polygon_vertices(polygon)
                            .is_some_and(|ids| self.oriented_chart_polygon_area(ids).is_some())
                    })
            })
        {
            return None;
        }
        let cuts = self.torus_cuts(&rim, &paths, ring)?;
        self.finish_variable_band(cycles, ring_index, ring, changed, cuts, trims)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shear_certifies_backwards_step_without_admitting_retracing() {
        let paths = [
            RimPath {
                advance: PI,
                height: [1.0, 1.0],
                kind: PathKind::Circle,
            },
            RimPath {
                advance: PI + 1e-10,
                height: [1.0, 1.0],
                kind: PathKind::Circle,
            },
            RimPath {
                advance: -1e-10,
                height: [1.0, 1.1],
                kind: PathKind::Lifted,
            },
        ];
        assert!(simple_cycle(&paths, 1e-8));
        let paths = [
            paths[0],
            paths[1],
            RimPath {
                advance: -1e-10,
                height: [1.1, 1.0],
                kind: PathKind::Lifted,
            },
            paths[2],
        ];
        assert!(!simple_cycle(&paths, 1e-8));
    }
}
