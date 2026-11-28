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
use monstertruck_geometry::proof_bezier_spans;
use std::{
    collections::BTreeMap,
    f64::consts::{PI, TAU},
};
pub(super) mod rational;
mod simplicity;

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

/// The actual oriented native edges, with their original trim ordinals.
pub(super) struct OrientedRim {
    pub(super) edges: Vec<Edge>,
    pub(super) sources: Vec<(usize, usize)>,
}

impl OrientedRim {
    pub(super) fn collect<'a>(
        cycle: &[OrientedTrim],
        mut lookup: impl FnMut(usize) -> Option<&'a [Edge]>,
    ) -> Option<Self> {
        let mut edges = Vec::new();
        let mut sources = Vec::new();
        for reference in cycle {
            let current = lookup(reference.trim)?;
            if reference.reversed {
                edges.extend(current.iter().rev().map(Edge::inverse));
                sources.extend(
                    (0..current.len())
                        .rev()
                        .map(|ordinal| (reference.trim, ordinal)),
                );
            } else {
                edges.extend(current.iter().cloned());
                sources.extend((0..current.len()).map(|ordinal| (reference.trim, ordinal)));
            }
        }
        if edges.is_empty()
            || edges
                .iter()
                .zip(edges.iter().cycle().skip(1))
                .any(|(a, b)| a.back() != b.front())
        {
            return None;
        }
        Some(Self { edges, sources })
    }

    fn reversed(&self, reverse: bool) -> Self {
        if reverse {
            Self {
                edges: self.edges.iter().rev().map(Edge::inverse).collect(),
                sources: self.sources.iter().rev().copied().collect(),
            }
        } else {
            Self {
                edges: self.edges.clone(),
                sources: self.sources.clone(),
            }
        }
    }
}

struct RimInputs {
    edges: Vec<Edge>,
    sources: Vec<(usize, usize)>,
    tensor_spans: Option<Vec<(RimPath, Point3)>>,
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
        || simplicity::simple_periodic_cycle(paths, resolution)
}

// A full meridional winding requires a ring torus. An azimuthal
// annulus may also occupy a regular positive-radius part of a horn or
// spindle torus, provided the entire meridian interval excludes poles.
pub(super) fn regular_radial_band(
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

    fn torus_path(
        &self,
        points: [Point3; 2],
        circle: Option<&Edge>,
        ring: Ring,
        side: f64,
    ) -> Option<RimPath> {
        let start = self.torus_chart(points[0])?;
        let end = self.torus_chart(points[1])?;
        let kind = if circle.is_some() {
            PathKind::Circle
        } else {
            PathKind::Lifted
        };
        let resolution = self.options.angular_resolution;
        let delta = match kind {
            PathKind::Lifted => {
                let delta = std::array::from_fn(|i| signed_delta(start[i], end[i]));
                if delta.iter().any(|d| d.abs() >= PI - resolution) {
                    return None;
                }
                delta
            }
            PathKind::Circle => self.circle_delta(circle?, start)?,
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

    // Lift provenance certifies each positive knot span is a complete
    // tensor-circle path. Extract all original spans without tolerant
    // knot grouping; a C0 spline is still one shared topological edge.
    fn torus_spans(
        &self,
        edge: &Edge,
        kind: PathKind,
        ring: Ring,
        side: f64,
    ) -> Option<Vec<(RimPath, Point3)>> {
        let curve = edge.oriented_curve();
        let front = edge.front().point();
        let back = edge.back().point();
        if matches!(kind, PathKind::Circle) {
            let (a, b) = curve.range_tuple();
            return Some(vec![(
                self.torus_path([front, back], Some(edge), ring, side)?,
                curve.evaluate(a.midpoint(b)),
            )]);
        }
        let spline = curve.try_into_homogeneous_bspline_curve()?;
        let spans = proof_bezier_spans(&spline)?;
        let tol = self.options.intersection_tolerance;
        let mut previous = front;
        let mut result = Vec::with_capacity(spans.len());
        let point = |&[x, y, z, w]: &[f64; 4]| Point3::new(x / w, y / w, z / w);
        for span in spans {
            if span.controls.len() != 5 {
                return None;
            }
            let a = point(span.controls.first()?);
            let b = point(span.controls.last()?);
            if (a - previous).norm() > tol {
                return None;
            }
            result.push((
                self.torus_path([a, b], None, ring, side)?,
                curve.evaluate(span.start.midpoint(span.end)),
            ));
            previous = b;
        }
        ((previous - back).norm() <= tol).then_some(result)
    }

    fn torus_cuts(
        &self,
        first_point: Point3,
        paths: &[RimPath],
        midpoints: &[Point3],
        ring: Ring,
    ) -> Option<[([f64; 2], Point3); 2]> {
        let resolution = self.options.angular_resolution;
        let direction = if ring.positive { -1.0 } else { 1.0 };
        let first = self.torus_chart(first_point)?[ring.axis];
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
        for (i, &point) in midpoints.iter().enumerate() {
            if paths[i].advance <= 2.0 * resolution {
                continue;
            }
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

    fn collect_rim(
        &self,
        cycle: &[OrientedTrim],
        changed: &BTreeMap<usize, Vec<Edge>>,
        ring: Ring,
        side: f64,
    ) -> Option<RimInputs> {
        let OrientedRim {
            edges: rim,
            sources: rim_sources,
        } = OrientedRim::collect(cycle, |id| changed.get(&id).map(Vec::as_slice))?;
        let mut spans = Vec::new();
        let mut tensor_only = true;
        for (edge, &(id, _)) in rim.iter().zip(&rim_sources) {
            let trim = self.graph.trims.get(id)?;
            let kind = if self.lifted_deviations.get(id).is_some_and(Option::is_some) {
                if changed.get(&id)?.len() != trim.source_edges.len() {
                    return None;
                }
                Some(PathKind::Lifted)
            } else if matches!(trim.geometry, SharedTrimGeometry::Circle(_)) {
                Some(PathKind::Circle)
            } else {
                None
            };
            if tensor_only {
                if let Some(parts) = kind.and_then(|kind| self.torus_spans(edge, kind, ring, side))
                {
                    spans.extend(parts);
                } else {
                    tensor_only = false;
                }
            }
        }
        Some(RimInputs {
            edges: rim,
            sources: rim_sources,
            tensor_spans: tensor_only.then_some(spans),
        })
    }

    fn tensor_rim_proof(&self, inputs: &RimInputs, ring: Ring) -> Option<rational::RimProof> {
        let (paths, midpoints): (Vec<_>, Vec<_>) =
            inputs.tensor_spans.as_ref()?.iter().copied().unzip();
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
        let cuts = self.torus_cuts(
            inputs.edges.first()?.front().point(),
            &paths,
            &midpoints,
            ring,
        )?;
        Some(rational::RimProof { extent, cuts })
    }

    pub(super) fn variable_torus_band(
        &self,
        cycles: [&[OrientedTrim]; 2],
        ring_index: usize,
        ring: Ring,
        trims: &mut [Result<Vec<Edge>, super::NativeTrimError>],
        obstacles: &mut Vec<super::NativeChartObstacle>,
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
        let inputs = self.collect_rim(cycles[1 - ring_index], &changed, ring, side)?;
        let rational::RimProof { extent, cuts } = match self.tensor_rim_proof(&inputs, ring) {
            Some(proof) => proof,
            None => match rational::prove(self, &inputs.edges, ring, side) {
                Ok(proof) => proof,
                Err(error) => {
                    obstacles.extend(chart_obstacles(error, &inputs.sources));
                    return None;
                }
            },
        };
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
            || !self.source_chart_positive()
        {
            obstacles.push(super::NativeChartObstacle::Coverage);
            return None;
        }
        let plan = self.finish_variable_band(cycles, ring_index, ring, changed, cuts, trims);
        if plan.is_none() {
            obstacles.push(super::NativeChartObstacle::SeamConstruction);
        }
        plan
    }
}

pub(super) fn chart_obstacles(
    error: rational::Rejection,
    sources: &[(usize, usize)],
) -> Vec<super::NativeChartObstacle> {
    use super::NativeChartObstacle as Obstacle;
    use rational::Rejection;
    match error {
        Rejection::Span(id) => vec![
            sources
                .get(id)
                .copied()
                .map_or(Obstacle::Invalid, Obstacle::Span),
        ],
        Rejection::Ring(id) => vec![
            sources
                .get(id)
                .copied()
                .map_or(Obstacle::Invalid, Obstacle::RingContact),
        ],
        Rejection::Pair(pairs) => {
            let mut result = pairs
                .into_iter()
                .filter_map(|[a, b]| Some(Obstacle::RimPair([*sources.get(a)?, *sources.get(b)?])))
                .collect::<Vec<_>>();
            if result.is_empty() {
                result.push(Obstacle::Invalid);
            }
            result
        }
        Rejection::Winding => vec![Obstacle::Winding],
        Rejection::Cuts => vec![Obstacle::CutSelection],
        Rejection::Invalid => vec![Obstacle::Invalid],
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
