//! Positive witnesses of missing source support on finite revolution faces.
//! The infinite carrier test cannot distinguish a correctly fitted surface
//! with a wrong trimming loop. Work in a regular lifted chart, using
//! interval enclosures of the actual rational boundary, not sampled chords.
//! An inconclusive chart never constitutes a missing-coverage witness.

use super::{BandContext, Frame, Interval, Piece, Ring, angle, lift_cycle, simplicity};
use crate::{Face, Point3, Wire, mesh_recovery::SurfaceCarrier};
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};

type ChartBox = [Interval; 2];

fn disjoint(a: ChartBox, b: ChartBox) -> bool {
    (0..2).any(|axis| a[axis].hi < b[axis].lo || b[axis].hi < a[axis].lo)
}

// If ||q-p|| <= budget, its axial/radial coordinates differ by at most
// that distance: (rho,z) is nonexpansive. The angular image of a disk of
// radius b around a vector of length r fits in ±asin(b/r). Bound asin(t)
// by t/sqrt(1-t²), avoiding another transcendental rounding assumption.
fn angular_ball(center: Interval, norm: Interval, budget: Interval) -> Option<Interval> {
    let relative = budget.div(norm)?;
    if relative.lo < 0.0 || relative.hi >= 1.0 {
        return None;
    }
    let denominator = Interval::point(1.0).sub(relative.square()).sqrt()?;
    let error = relative.div(denominator)?.hi;
    Some(center.add(Interval::range(-error, error)))
}

fn chart_budget(frame: Frame, budget: f64) -> Option<Interval> {
    // The local basis is orthonormal up to rounding. Include its measured
    // operator norm upper bound and arithmetic roundoff, not a geometric
    // tolerance expansion. For orthonormal axes, the Gershgorin bound is 1.
    let norm2 = frame
        .basis
        .iter()
        .map(|a| {
            frame
                .basis
                .iter()
                .map(|b| {
                    (0..3)
                        .map(|i| Interval::point(a[i]).mul(Interval::point(b[i])))
                        .fold(Interval::point(0.0), Interval::add)
                        .magnitude()
                })
                .map(Interval::point)
                .fold(Interval::point(0.0), Interval::add)
                .hi
        })
        .fold(0.0_f64, f64::max);
    Some(Interval::point(budget).mul(Interval::point(norm2).sqrt()?))
}

fn source_ball(frame: Frame, p: Point3, budget: Interval) -> Option<ChartBox> {
    let [x, y, z, _] = frame.local([p.x, p.y, p.z, 1.0].map(Interval::point));
    let rho = x.square().add(y.square()).sqrt()?;
    let u = angular_ball(angle(x, y)?, rho, budget)?;
    let v = if let Some(major) = frame.major {
        let radial = rho.sub(Interval::point(major));
        let norm = radial.square().add(z.square()).sqrt()?;
        angular_ball(angle(radial, z)?, norm, budget)?
    } else {
        z.add(Interval::range(-budget.hi, budget.hi))
    };
    Some(frame.orient([u, v]))
}

struct Domain {
    pieces: Vec<Piece>,
    bounds: ChartBox,
}

impl Domain {
    fn new(context: &BandContext<'_>, frame: Frame, wire: &Wire) -> Option<Self> {
        let edges = wire.iter().cloned().collect::<Vec<_>>();
        if edges
            .iter()
            .zip(edges.iter().cycle().skip(1))
            .any(|(a, b)| a.back() != b.front())
        {
            return None;
        }
        let first = frame.endpoint(frame.local({
            let p = edges.first()?.front().point();
            [p.x, p.y, p.z, 1.0].map(Interval::point)
        }))?;
        let (pieces, end) = lift_cycle(frame, &edges, first).ok()?;
        // A noncontractible or singular wire needs a different chart.
        if (0..2).any(|i| end[i].sub(first[i]).magnitude() > frame.tolerance(i)) {
            return None;
        }
        // Tolerance-close endpoints do not by themselves establish a
        // continuous chart cycle. Require compatible interval endpoints
        // at every polynomial span junction, including the closing seam.
        if pieces
            .iter()
            .zip(pieces.iter().cycle().skip(1))
            .any(|(a, b)| (0..2).any(|axis| a.end[axis].intersection(b.start[axis]).is_none()))
        {
            return None;
        }
        let bounds = pieces
            .iter()
            .map(|p| p.bounds)
            .reduce(|a, b| std::array::from_fn(|i| a[i].union(b[i])))?;
        if bounds[0].width() >= TAU || frame.periodic_height() && bounds[1].width() >= TAU {
            return None;
        }
        // The entire compact outer chart must lie on the positive radial
        // sheet, not merely its boundary. Otherwise the canonical angular
        // map is not injective on the purported bounded face domain.
        let regular = match context.carrier {
            SurfaceCarrier::CircularRevolution {
                radial_center,
                profile_radius,
                ..
            } => {
                let v = bounds[1];
                let minimum = if (PI - v.lo).rem_euclid(TAU) <= v.width() {
                    -1.0
                } else {
                    v.lo.cos().min(v.hi.cos())
                };
                Interval::point(minimum)
                    .sub(Interval::point(128.0 * f64::EPSILON))
                    .scale(profile_radius)
                    .add(Interval::point(radial_center))
                    .lo
                    > context.options.intersection_tolerance
            }
            SurfaceCarrier::LinearRevolution {
                radius_at_origin,
                slope,
                ..
            } => {
                bounds[1]
                    .scale(slope)
                    .add(Interval::point(radius_at_origin))
                    .lo
                    > context.options.intersection_tolerance
            }
            SurfaceCarrier::Plane { .. } => false,
        };
        if !regular || simplicity::cycle(frame, &pieces, false).is_err() {
            return None;
        }
        // Certify the oriented compact side of the loop, rather than
        // mistaking a complement on a periodic surface for its interior.
        // On each normalized Bezier span, the interval integrand encloses
        // x dy - y dx for the complete curved boundary.
        let area = pieces
            .iter()
            .map(|p| {
                let x = p.bounds[0].sub(first[0]);
                let y = p.bounds[1].sub(first[1]);
                x.mul(p.rates[1]).sub(y.mul(p.rates[0]))
            })
            .fold(Interval::point(0.0), Interval::add);
        if area.scale(if context.outward { 1.0 } else { -1.0 }).lo <= 0.0 {
            return None;
        }
        Some(Self { pieces, bounds })
    }

    fn excludes(&self, frame: Frame, query: ChartBox) -> bool {
        if query[0].width() >= TAU || frame.periodic_height() && query[1].width() >= TAU {
            return false;
        }
        // A compact lifted chart may straddle either periodic seam. Test
        // every potentially overlapping copy, not just the nearest center.
        let query: ChartBox =
            std::array::from_fn(|i| frame.unwrap(i, query[i], self.bounds[i].mid()));
        let heights: &[f64] = if frame.periodic_height() {
            &[-TAU, 0.0, TAU]
        } else {
            &[0.0]
        };
        [-TAU, 0.0, TAU].into_iter().all(|u| {
            heights.iter().all(|&v| {
                self.excludes_image([
                    query[0].add(Interval::point(u)),
                    query[1].add(Interval::point(v)),
                ])
            })
        })
    }

    fn excludes_image(&self, query: ChartBox) -> bool {
        if disjoint(self.bounds, query) {
            return true;
        }
        // Each rational span and its endpoint chord are contained in the
        // same convex chart box. If the entire query box is disjoint from
        // every span, chord homotopy preserves winding around any of its
        // points, even for nonmonotone trims. A ray through an endpoint is
        // inconclusive; try the orthogonal ray as well.
        if self
            .pieces
            .iter()
            .any(|piece| !disjoint(piece.bounds, query))
        {
            return false;
        }
        let center = query.map(Interval::mid);
        (0..2).any(|axis| {
            self.pieces.iter().try_fold(0_i32, |sum, piece| {
                Some(sum + ray_crossing(piece, center, axis)?)
            }) == Some(0)
        })
    }
}

// Use interval signs, not rounded chart endpoint midpoints. The ray is
// positive in axis x; swapping the coordinate roles only changes the
// winding sign, which does not affect a zero-winding witness.
fn ray_crossing(piece: &Piece, center: [f64; 2], x: usize) -> Option<i32> {
    let y = 1 - x;
    let [qx, qy] = [center[x], center[y]].map(Interval::point);
    let a = piece.start[y].sub(qy);
    let b = piece.end[y].sub(qy);
    if a.hi < 0.0 && b.hi < 0.0 || a.lo > 0.0 && b.lo > 0.0 || piece.bounds[x].hi < qx.lo {
        return Some(0);
    }
    let direction = if a.hi < 0.0 && b.lo > 0.0 {
        1
    } else if a.lo > 0.0 && b.hi < 0.0 {
        -1
    } else {
        return None;
    };
    let dx = piece.end[x].sub(piece.start[x]);
    let dy = piece.end[y].sub(piece.start[y]);
    let determinant = dx
        .mul(qy.sub(piece.start[y]))
        .sub(dy.mul(qx.sub(piece.start[x])))
        .scale(f64::from(direction));
    if determinant.lo > 0.0 {
        Some(direction)
    } else if determinant.hi < 0.0 {
        Some(0)
    } else {
        None
    }
}

pub(in super::super::super::super) fn missing_source_vertex(
    context: &BandContext<'_>,
    faces: &[Face],
    budget: f64,
) -> Option<usize> {
    let frame = Frame::new(
        context,
        Ring {
            axis: 0,
            fixed: 0.0,
            positive: false,
        },
        1.0,
    )?;
    let budget = chart_budget(frame, budget)?;
    // Only the outer loop is needed for a positive exclusion witness.
    // Ignoring holes may miss exclusions, but cannot create one.
    let domains = faces
        .iter()
        .map(|face| Domain::new(context, frame, face.boundaries().first()?))
        .collect::<Option<Vec<_>>>()?;
    if domains.is_empty() {
        return None;
    }
    let mut vertices = BTreeSet::new();
    for &support in context.supports {
        for &polygon in &context
            .evidence
            .planar_supports
            .get(support)?
            .source_face_indices
        {
            vertices.extend(
                context
                    .evidence
                    .source_polygon_vertices(polygon)?
                    .iter()
                    .copied(),
            );
        }
    }
    vertices.into_iter().find(|&id| {
        context
            .evidence
            .positions
            .get(id)
            .and_then(|&point| source_ball(frame, point, budget))
            .is_some_and(|query| domains.iter().all(|domain| domain.excludes(frame, query)))
    })
}

#[cfg(test)]
mod tests {
    use super::super::Controls;
    use super::*;

    fn frame() -> Frame {
        Frame {
            center: Point3::new(0.0, 0.0, 0.0),
            basis: [
                crate::Vector3::unit_x(),
                crate::Vector3::unit_y(),
                crate::Vector3::unit_z(),
            ],
            major: None,
            varying: 0,
            direction: 1.0,
            side: 1.0,
            resolution: 1e-8,
            height_resolution: 1e-8,
        }
    }

    fn polygon(points: &[[f64; 2]]) -> Option<Domain> {
        let pieces = points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .enumerate()
            .map(|(edge, (a, b))| {
                let controls = Controls {
                    values: vec![[Interval::point(0.0); 4]; 2],
                    parameter: [0.0, 1.0],
                };
                Piece {
                    controls,
                    start: a.map(Interval::point),
                    end: b.map(Interval::point),
                    bounds: std::array::from_fn(|i| Interval::range(a[i], b[i])),
                    rates: std::array::from_fn(|i| Interval::point(b[i] - a[i])),
                    edge,
                }
            })
            .collect::<Vec<_>>();
        let bounds = pieces
            .iter()
            .map(|p| p.bounds)
            .reduce(|a, b| std::array::from_fn(|i| a[i].union(b[i])))?;
        Some(Domain { pieces, bounds })
    }

    #[test]
    fn local_exclusion_inside_global_parameter_bounds() -> Result<(), Box<dyn std::error::Error>> {
        // A notch excludes the interior witness even though both of its
        // coordinates lie strictly inside the overall outer loop bounds.
        let domain = polygon(&[
            [0., 0.],
            [2., 0.],
            [2., 2.],
            [1.2, 2.],
            [1.2, 0.8],
            [0.8, 0.8],
            [0.8, 2.],
            [0., 2.],
        ])
        .ok_or("nonempty fixture")?;
        let query = |u, v| {
            [
                Interval::range(u - 0.05, u + 0.05),
                Interval::range(v - 0.05, v + 0.05),
            ]
        };
        assert!(domain.excludes(frame(), query(1., 1.5)));
        assert!(!domain.excludes(frame(), query(0.5, 1.5)));
        assert!(!domain.excludes(frame(), query(0.8, 1.5)));
        assert!(!domain.excludes(frame(), query(0.5 - TAU, 1.5)));
        assert!(!domain.excludes(frame(), [Interval::range(-TAU, TAU), Interval::point(1.5)]));
        Ok(())
    }

    #[test]
    fn source_ball_bounds_angular_and_axial_neighborhood() -> Result<(), Box<dyn std::error::Error>>
    {
        let frame = frame();
        let point = Point3::new(1.0, 0.0, 2.0);
        let bounds = source_ball(
            frame,
            point,
            chart_budget(frame, 0.1).ok_or("chart budget")?,
        )
        .ok_or("regular chart")?;
        assert!(bounds[0].contains(0.1));
        assert!(bounds[1].contains(2.1));
        assert!(
            source_ball(
                frame,
                point,
                chart_budget(frame, 1.1).ok_or("chart budget")?
            )
            .is_none()
        );
        Ok(())
    }
}
