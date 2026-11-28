//! Annular chart decomposition with contractible interior boundaries.
use super::{
    BandContext, BandPlan, NativeChartObstacle, NativeTrimError, RegionTrimLoops, SurfaceCarrier,
    torus_rim,
};
use crate::Edge;
use std::f64::consts::{PI, TAU};

type Bounds = [[f64; 2]; 2];

fn disjoint(a: Bounds, b: Bounds, tol: f64) -> bool {
    if a[1][1] + tol < b[1][0] || b[1][1] + tol < a[1][0] {
        return true;
    }
    let midpoint = |r: [f64; 2]| r[0].midpoint(r[1]);
    let shift = ((midpoint(a[0]) - midpoint(b[0])) / TAU).round() * TAU;
    [-TAU, 0.0, TAU].into_iter().all(|image| {
        let offset = shift + image;
        a[0][1] + tol < b[0][0] + offset || b[0][1] + offset + tol < a[0][0]
    })
}

fn clear(angle: f64, bound: [f64; 2], tol: f64) -> bool {
    let lifted = ((bound[0].midpoint(bound[1]) - angle) / TAU)
        .round()
        .mul_add(TAU, angle);
    lifted < bound[0] - tol || lifted > bound[1] + tol
}

// A seam is a full constant-azimuth/meridian connector. Exclude its
// intersection with every certified hole box. Prefer antipodal cuts,
// then choose the most separated clear gap midpoints if necessary.
fn cuts(holes: &[(usize, Bounds)], tol: f64) -> Option<[f64; 2]> {
    let mut events = vec![0.0, TAU];
    for (_, b) in holes {
        for angle in [b[0][0] - tol, b[0][1] + tol] {
            events.push(angle.rem_euclid(TAU));
            events.push((angle + PI).rem_euclid(TAU));
        }
    }
    events.sort_by(f64::total_cmp);
    let admissible = |angle| holes.iter().all(|(_, b)| clear(angle, b[0], tol));
    let candidates = events
        .windows(2)
        .filter(|p| p[1] - p[0] > 3.0 * tol)
        .flat_map(|p| {
            [
                p[0].midpoint(p[1]),
                (p[1] - p[0]).mul_add(1.0 / 3.0, p[0]),
                (p[1] - p[0]).mul_add(2.0 / 3.0, p[0]),
            ]
        })
        .filter(|&a| admissible(a))
        .collect::<Vec<_>>();
    for &a in &candidates {
        let b = (a + PI).rem_euclid(TAU);
        if admissible(b) {
            return Some([a, b]);
        }
    }
    candidates
        .iter()
        .flat_map(|&a| candidates.iter().map(move |&b| [a, b]))
        .max_by(|a, b| angular_distance(*a).total_cmp(&angular_distance(*b)))
        .filter(|&pair| angular_distance(pair) > tol)
}

fn angular_distance([a, b]: [f64; 2]) -> f64 { ((a - b + PI).rem_euclid(TAU) - PI).abs() }

impl BandContext<'_> {
    pub(super) fn holed_band(
        &self,
        boundary: &RegionTrimLoops,
        trims: &mut [Result<Vec<Edge>, NativeTrimError>],
        obstacles: &mut Vec<NativeChartObstacle>,
    ) -> Option<BandPlan> {
        let SurfaceCarrier::CircularRevolution {
            radial_center,
            profile_radius,
            ..
        } = self.carrier
        else {
            return None;
        };
        let rings = boundary
            .loops
            .iter()
            .enumerate()
            .filter_map(|(i, cycle)| self.ring(cycle).map(|r| (i, r)))
            .collect::<Vec<_>>();
        let [(first, a), (second, b)] = rings.as_slice() else {
            return None;
        };
        let rings = [*a, *b];
        let (first_lower, lower, span) = self.ring_domain(rings)?;
        let axis = a.axis;
        let tol = self.options.angular_resolution;
        if !torus_rim::regular_radial_band(
            radial_center,
            profile_radius,
            axis,
            lower,
            span,
            self.options.intersection_tolerance,
        ) || !self.covers_band(1 - axis, lower, span)
            || !self.source_chart_positive()
        {
            return None;
        }
        let mut holes: Vec<(usize, Bounds)> = Vec::new();
        let mut rejected = false;
        for (i, cycle) in boundary.loops.iter().enumerate() {
            if i == *first || i == *second {
                continue;
            }
            let mut rim = Vec::new();
            let mut sources = Vec::new();
            for reference in cycle {
                let edges = trims.get(reference.trim)?.as_ref().ok()?;
                if reference.reversed {
                    sources.extend(
                        (0..edges.len())
                            .rev()
                            .map(|ordinal| (reference.trim, ordinal)),
                    );
                    rim.extend(edges.iter().rev().map(Edge::inverse));
                } else {
                    sources.extend((0..edges.len()).map(|ordinal| (reference.trim, ordinal)));
                    rim.extend(edges.iter().cloned());
                }
            }
            let bounds = match torus_rim::rational::hole_bounds(self, &rim, axis, lower, span) {
                Ok(bounds) => bounds,
                Err(error) => {
                    rejected = true;
                    obstacles.extend(torus_rim::chart_obstacles(error, &sources));
                    continue;
                }
            };
            // Pairwise disjoint boxes certify both disjointness and lack
            // of nesting; otherwise no hole admission is asserted.
            if holes
                .iter()
                .any(|(_, other)| !disjoint(bounds, *other, tol))
            {
                return None;
            }
            holes.push((i, bounds));
        }
        if rejected {
            return None;
        }
        let angles = cuts(&holes, tol).or_else(|| {
            obstacles.push(NativeChartObstacle::CutSelection);
            None
        })?;
        let sign = if a.positive { 1.0 } else { -1.0 };
        let sweep = (sign * (angles[1] - angles[0])).rem_euclid(TAU);
        let mut assigned = [Vec::new(), Vec::new()];
        for (i, bound) in holes {
            let phase = (sign * (bound[0][0].midpoint(bound[0][1]) - angles[0])).rem_euclid(TAU);
            let index = usize::from(phase > sweep);
            assigned[index].push(i);
        }
        let mut plan = self.parallel_band(
            [&boundary.loops[*first], &boundary.loops[*second]],
            rings,
            angles,
            if first_lower { span } else { -span },
            trims,
        )?;
        plan.ring_wires = [*first, *second];
        plan.holes = assigned;
        Some(plan)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cuts_clear_wrapped_holes_and_do_not_require_antipodal_gaps() {
        let tolerance = 1.0e-8;
        let holes = [
            (2, [[-0.5, 0.5], [0.1, 0.2]]),
            (3, [[1.0, 5.5], [0.3, 0.4]]),
        ];
        let angles = cuts(&holes, tolerance).expect("separate clear seam angles");
        assert!(angular_distance(angles) > tolerance);
        assert!(
            angles
                .iter()
                .all(|&angle| holes.iter().all(|(_, b)| clear(angle, b[0], tolerance)))
        );
        assert!(!disjoint(
            [[-0.1, 0.3], [0.1, 0.2]],
            [[TAU - 0.2, TAU + 0.1], [0.15, 0.25]],
            tolerance
        ));
    }
}
