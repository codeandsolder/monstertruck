//! A regular linear-revolution band with one exact isoparametric rim and
//! one piecewise varying-height rim. Prove the latter is an azimuth-monotone
//! cycle, including vertical step edges, before inserting shared chart cuts.

use super::{
    BandContext, BandPlan, Ring, commit_rings, coordinate_difference, copy_rings, cut_ring,
    world_point,
};
use crate::mesh_recovery::{OrientedTrim, SurfaceCarrier};
use crate::{BoundedCurve, Curve, Edge, InnerSpace, ParametricCurve, Point3, Vector3, Vertex};
use monstertruck_geometry::{
    prelude::TryIntoHomogeneousBsplineCurve, proof_bezier_spans, subdivide_homogeneous_bezier,
};
use std::f64::consts::{PI, TAU};

struct RimEdge {
    advance: f64,
    axial: [f64; 2],
}

struct RimChart {
    origin: Point3,
    axis: Vector3,
    sign: f64,
    tolerance: f64,
    angular_resolution: f64,
}

impl RimChart {
    fn radial(&self, point: Point3) -> Vector3 {
        let delta = point - self.origin;
        delta - self.axis * delta.dot(self.axis)
    }

    fn angle(&self, a: Vector3, b: Vector3) -> f64 {
        self.sign * self.axis.dot(a.cross(b)).atan2(a.dot(b))
    }

    // A positive rational combination of ordered radial control vectors
    // in one open half-plane has monotone azimuth. All pairwise determinants
    // then have the same sign; the rational derivative has that sign too.
    // Recursively split wide spans rather than interpreting an ambiguous
    // half-turn as the desired branch.
    fn span(&self, controls: &[[f64; 4]], depth: usize) -> Option<RimEdge> {
        let mut axial = [f64::INFINITY, f64::NEG_INFINITY];
        let mut directions = Vec::with_capacity(controls.len());
        let mut scale = 1.0_f64;
        let mut weights = [f64::INFINITY, 0.0_f64];
        for &[x, y, z, w] in controls {
            if !w.is_finite() || w <= 0.0 {
                return None;
            }
            let point = Point3::new(x / w, y / w, z / w);
            let radial = self.radial(point);
            let radius = radial.magnitude();
            if !radius.is_finite() || radius <= self.tolerance {
                return None;
            }
            let height = (point - self.origin).dot(self.axis);
            if !height.is_finite() {
                return None;
            }
            scale = scale
                .max((point - self.origin).magnitude())
                .max(point.x.abs())
                .max(point.y.abs())
                .max(point.z.abs());
            weights[0] = weights[0].min(w);
            weights[1] = weights[1].max(w);
            directions.push((radial / radius, radius));
            axial[0] = axial[0].min(height);
            axial[1] = axial[1].max(height);
        }
        let degree = crate::usize_to_f64(controls.len());
        let margin = 128.0 * f64::EPSILON * scale * (weights[1] / weights[0]) * degree * degree;
        if !margin.is_finite() {
            return None;
        }
        let first = directions.first()?.0;
        let last = directions.last()?.0;
        let advance = self.angle(first, last);
        let mut previous = first;
        let mut monotone =
            advance > self.angular_resolution && advance < PI - self.angular_resolution;
        for &(direction, radius) in directions.iter().skip(1) {
            let offset = self.angle(first, direction);
            let step = self.angle(previous, direction);
            // Require a strict sign outside the arithmetic uncertainty
            // unless the normalized direction is represented identically.
            monotone &= offset >= 0.0
                && offset <= advance
                && (previous == direction || step > margin / radius);
            previous = direction;
        }
        if monotone {
            return Some(RimEdge {
                advance,
                axial: [axial[0] - margin, axial[1] + margin],
            });
        }
        if depth == 0 {
            return None;
        }
        let (a, b) = subdivide_homogeneous_bezier(controls)?;
        let a = self.span(&a, depth - 1)?;
        let b = self.span(&b, depth - 1)?;
        Some(RimEdge {
            advance: a.advance + b.advance,
            axial: [a.axial[0].min(b.axial[0]), a.axial[1].max(b.axial[1])],
        })
    }

    fn edge(&self, edge: &Edge) -> Option<RimEdge> {
        let curve = edge.oriented_curve();
        if matches!(curve, Curve::Line(_)) {
            let a = edge.front().point();
            let b = edge.back().point();
            let ar = self.radial(a);
            let br = self.radial(b);
            // A line on a regular linear revolution is a generator if both
            // endpoints lie on the same radial ray. Its whole axial range
            // is given by those endpoints.
            if ar.dot(br) <= 0.0
                || ar.cross(br).magnitude() > self.tolerance * ar.magnitude().max(br.magnitude())
            {
                return None;
            }
            let z0 = (a - self.origin).dot(self.axis);
            let z1 = (b - self.origin).dot(self.axis);
            if (z1 - z0).abs() <= self.tolerance {
                return None;
            }
            return Some(RimEdge {
                advance: 0.0,
                axial: [z0.min(z1), z0.max(z1)],
            });
        }
        let curve = curve.try_into_homogeneous_bspline_curve()?;
        let mut result = RimEdge {
            advance: 0.0,
            axial: [f64::INFINITY, f64::NEG_INFINITY],
        };
        for span in proof_bezier_spans(&curve)? {
            let span = self.span(&span.controls, 8)?;
            result.advance += span.advance;
            result.axial[0] = result.axial[0].min(span.axial[0]);
            result.axial[1] = result.axial[1].max(span.axial[1]);
        }
        (result.advance > 0.0).then_some(result)
    }
}

fn nonretracing_rim(rim: &[Edge], proofs: &[RimEdge], axis: Vector3) -> bool {
    // Consecutive generators at one azimuth must not retrace the rim.
    (0..rim.len()).all(|i| {
        let j = (i + 1) % rim.len();
        proofs[i].advance != 0.0
            || proofs[j].advance != 0.0
            || (rim[i].back().point() - rim[i].front().point()).dot(axis)
                * (rim[j].back().point() - rim[j].front().point()).dot(axis)
                > 0.0
    })
}

fn rim_extent(
    rim: &[Edge],
    proofs: &[RimEdge],
    axis: Vector3,
    angular_resolution: f64,
) -> Option<[f64; 2]> {
    let sweep = proofs.iter().map(|proof| proof.advance).sum::<f64>();
    if (sweep - TAU).abs() > angular_resolution || !nonretracing_rim(rim, proofs, axis) {
        return None;
    }
    Some(
        proofs
            .iter()
            .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], proof| {
                [lo.min(proof.axial[0]), hi.max(proof.axial[1])]
            }),
    )
}

impl BandContext<'_> {
    pub(super) fn variable_band(
        &self,
        cycles: [&[OrientedTrim]; 2],
        ring_index: usize,
        ring: Ring,
        trims: &mut [Result<Vec<Edge>, super::NativeTrimError>],
    ) -> Option<BandPlan> {
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } = self.carrier
        else {
            return None;
        };
        if ring.axis != 0 {
            return None;
        }
        let other = 1 - ring_index;
        let mut changed = copy_rings(cycles[0].iter().chain(cycles[1]), trims)?;
        let mut rim = Vec::new();
        for reference in cycles[other] {
            let edges = changed.get(&reference.trim)?;
            if reference.reversed {
                rim.extend(edges.iter().rev().map(Edge::inverse));
            } else {
                rim.extend(edges.iter().cloned());
            }
        }
        let chart = RimChart {
            origin,
            axis,
            sign: if ring.positive { -1.0 } else { 1.0 },
            tolerance: self.options.intersection_tolerance,
            angular_resolution: self.options.angular_resolution,
        };
        let proofs = rim
            .iter()
            .map(|edge| chart.edge(edge))
            .collect::<Option<Vec<_>>>()?;
        let [lo, hi] = rim_extent(&rim, &proofs, axis, self.options.angular_resolution)?;
        let ring_lower = ring.positive == self.outward;
        let tol = self.options.intersection_tolerance;
        if (ring_lower && lo <= ring.fixed + tol) || (!ring_lower && hi >= ring.fixed - tol) {
            return None;
        }
        let lower = lo.min(ring.fixed);
        let upper = hi.max(ring.fixed);
        if slope
            .mul_add(lower, radius_at_origin)
            .min(slope.mul_add(upper, radius_at_origin))
            <= tol
            || !self.covers_band(1, lower, upper - lower)
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
        let [first, second] = self.rim_cuts(&rim, &proofs)?;
        let mut vertices = [Vec::<Vertex>::new(), Vec::new()];
        let mut seams = Vec::new();
        for (uv, rim_point) in [first, second] {
            let point = world_point(self.carrier, [uv[0], ring.fixed])?;
            let ring_vertex = cut_ring(cycles[ring_index], point, &mut changed, tol)?;
            let rim_vertex = cut_ring(cycles[other], rim_point, &mut changed, tol)?;
            vertices[ring_index].push(ring_vertex);
            vertices[other].push(rim_vertex);
            let (start, end) = if ring_index == 0 {
                (ring.fixed, uv[1])
            } else {
                (uv[1], ring.fixed)
            };
            seams.push(self.seam(
                0,
                uv[0],
                start,
                end - start,
                vertices[0].last()?,
                vertices[1].last()?,
            )?);
        }
        let vertices = vertices.map(|v| v.try_into().ok());
        let [Some(a), Some(b)] = vertices else {
            return None;
        };
        let seams = seams.try_into().ok()?;
        commit_rings(changed, trims);
        Some(BandPlan {
            vertices: [a, b],
            seams,
        })
    }
    fn rim_cuts(&self, rim: &[Edge], proofs: &[RimEdge]) -> Option<[([f64; 2], Point3); 2]> {
        // Prefer existing native junctions with advancing azimuth on both
        // sides. If the cycle alternates arcs and axial steps, subdivide
        // the actual analytic arcs at interior parameters instead.
        let project = |point| {
            super::chart_point(self.carrier, point, self.options.angular_resolution)
                .map(|uv| (uv, point))
        };
        let mut candidates = rim
            .iter()
            .enumerate()
            .filter_map(|(i, edge)| {
                let previous = (i + rim.len() - 1) % rim.len();
                (proofs[i].advance > 0.0 && proofs[previous].advance > 0.0)
                    .then(|| project(edge.front().point()))
                    .flatten()
            })
            .collect::<Vec<_>>();
        if candidates.len() < 2 {
            candidates.extend(rim.iter().zip(proofs).filter_map(|(edge, proof)| {
                if proof.advance <= 0.0 {
                    return None;
                }
                let curve = edge.oriented_curve();
                let (a, b) = curve.range_tuple();
                let point = curve.evaluate(a.midpoint(b));
                let clearance = 4.0 * self.options.intersection_tolerance;
                ((point - edge.front().point()).magnitude() > clearance
                    && (point - edge.back().point()).magnitude() > clearance)
                    .then(|| project(point))
                    .flatten()
            }));
        }
        let first = *candidates.first()?;
        let distance = |a| coordinate_difference(self.carrier, 0, a, first.0[0]).abs();
        let second = *candidates
            .iter()
            .max_by(|a, b| distance(a.0[0]).total_cmp(&distance(b.0[0])))?;
        (distance(second.0[0]) > self.options.angular_resolution).then_some([first, second])
    }
}
