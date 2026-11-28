//! Positive witnesses of missing source support on finite revolution faces.
//! The infinite carrier test cannot distinguish a correctly fitted surface
//! with a wrong trimming loop. Work in a regular lifted chart, using
//! interval enclosures of the actual rational boundary, not sampled chords.
//! An inconclusive chart never constitutes a missing-coverage witness.

use super::{BandContext, Chart, Frame, Interval, Piece, Ring, angle, lift_cycle, simplicity};
use crate::{Face, InnerSpace, Point3, Vector3, Wire, mesh_recovery::SurfaceCarrier};
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

fn source_ball(frame: Frame, point: Point3, budget: Interval) -> Option<ChartBox> {
    let [x, y, z, _] = frame.local([point.x, point.y, point.z, 1.0].map(Interval::point));
    if matches!(frame.chart, Chart::Orthographic { .. }) {
        let extent = Interval::range(-budget.hi, budget.hi);
        return Some(frame.orient([x.add(extent), y.add(extent)]));
    }
    let rho = x.square().add(y.square()).sqrt()?;
    let azimuth = angular_ball(angle(x, y)?, rho, budget)?;
    let meridian = if let Some(major) = frame.major {
        let radial = rho.sub(Interval::point(major));
        let norm = radial.square().add(z.square()).sqrt()?;
        angular_ball(angle(radial, z)?, norm, budget)?
    } else {
        z.add(Interval::range(-budget.hi, budget.hi))
    };
    Some(frame.orient([azimuth, meridian]))
}

fn regular_revolution(context: &BandContext<'_>, bounds: ChartBox) -> bool {
    match context.carrier {
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
    }
}

// A wide span enclosure can obscure the orientation of a narrow but
// valid face. Subdivide only until the signed area has a definite sign.
fn oriented_pieces(
    frame: Frame,
    mut pieces: Vec<Piece>,
    origin: ChartBox,
    positive: bool,
) -> Option<Vec<Piece>> {
    let sign = if positive { 1.0 } else { -1.0 };
    loop {
        let area = pieces
            .iter()
            .map(|p| {
                let x = p.bounds[0].sub(origin[0]);
                let y = p.bounds[1].sub(origin[1]);
                x.mul(p.rates[1]).sub(y.mul(p.rates[0]))
            })
            .fold(Interval::point(0.0), Interval::add)
            .scale(sign);
        if area.lo > 0.0 {
            return Some(pieces);
        }
        if area.hi <= 0.0 || pieces.len() > super::MAX_PIECES / 2 {
            return None;
        }
        pieces = pieces
            .iter()
            .map(|piece| piece.split(frame))
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect();
    }
}

struct Domain {
    pieces: Vec<Piece>,
    bounds: ChartBox,
}

impl Domain {
    fn new(context: &BandContext<'_>, frame: Frame, wire: &Wire, positive: bool) -> Option<Self> {
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
        if (0..2).any(|axis| frame.periodic(axis) && bounds[axis].width() >= TAU) {
            return None;
        }
        // The entire compact outer chart must lie on the positive radial
        // sheet, not merely its boundary. Otherwise the canonical angular
        // map is not injective on the purported bounded face domain.
        let regular = match frame.chart {
            Chart::Orthographic { .. } => true,
            Chart::Revolution => regular_revolution(context, bounds),
        };
        if !regular || simplicity::cycle(frame, &pieces, false).is_err() {
            return None;
        }
        // Certify the oriented compact side of the loop, rather than
        // mistaking a complement on a periodic surface for its interior.
        // On each normalized Bezier span, the interval integrand encloses
        // x dy - y dx for the complete curved boundary.
        let pieces = oriented_pieces(frame, pieces, first, positive)?;
        Some(Self { pieces, bounds })
    }

    fn excludes(&self, frame: Frame, query: ChartBox) -> bool {
        if (0..2).any(|axis| frame.periodic(axis) && query[axis].width() >= TAU) {
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
        let azimuths: &[f64] = if frame.periodic(0) {
            &[-TAU, 0.0, TAU]
        } else {
            &[0.0]
        };
        azimuths.iter().all(|&u| {
            heights.iter().all(|&v| {
                self.excludes_image(
                    frame,
                    [
                        query[0].add(Interval::point(u)),
                        query[1].add(Interval::point(v)),
                    ],
                )
            })
        })
    }

    fn excludes_image(&self, frame: Frame, query: ChartBox) -> bool {
        if disjoint(self.bounds, query) {
            return true;
        }
        // Each rational span and its endpoint chord are contained in the
        // same convex chart box. If the entire query box is disjoint from
        // every span, chord homotopy preserves winding around any of its
        // points, even for nonmonotone trims. A ray through an endpoint is
        // inconclusive; try the orthogonal ray as well.
        let center = query.map(Interval::mid);
        (0..2).any(|axis| {
            let mut remaining = super::MAX_PIECES;
            self.pieces.iter().try_fold(0_i32, |sum, piece| {
                Some(
                    sum + separated_crossing(
                        frame,
                        piece,
                        query,
                        center,
                        axis,
                        super::MAX_DEPTH,
                        &mut remaining,
                    )?,
                )
            }) == Some(0)
        })
    }
}

// Recursively separate the actual rational span from the whole source
// ball's chart enclosure, not merely from its midpoint. Once separated,
// replacing the span by its endpoint chord preserves winding.
fn separated_crossing(
    frame: Frame,
    piece: &Piece,
    query: ChartBox,
    center: [f64; 2],
    axis: usize,
    depth: usize,
    remaining: &mut usize,
) -> Option<i32> {
    *remaining = remaining.checked_sub(1)?;
    if disjoint(piece.bounds, query) {
        return ray_crossing(piece, center, axis);
    }
    let inside =
        |point: ChartBox| (0..2).all(|i| query[i].lo <= point[i].lo && point[i].hi <= query[i].hi);
    if depth == 0 || inside(piece.start) || inside(piece.end) {
        return None;
    }
    let [first, last] = piece.split(frame)?;
    Some(
        separated_crossing(frame, &first, query, center, axis, depth - 1, remaining)?
            + separated_crossing(frame, &last, query, center, axis, depth - 1, remaining)?,
    )
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

// A single-sheet cone is a graph over its transverse plane, including
// the apex. A sphere restricted to an open hemisphere is likewise a
// graph; the sign of the projected boundary selects the compact side.
// These charts do not introduce angular singularities or artificial cuts.
fn projection_frame(
    context: &BandContext<'_>,
    origin: Point3,
    axis: Vector3,
    hemisphere: bool,
) -> Frame {
    let (u, v) = super::super::super::normal_basis(axis);
    Frame {
        chart: Chart::Orthographic { hemisphere },
        center: origin,
        basis: [u, v, axis],
        major: None,
        varying: 0,
        direction: 1.0,
        side: 1.0,
        resolution: context.options.intersection_tolerance,
        height_resolution: context.options.intersection_tolerance,
    }
}

struct FaceDomain {
    frame: Frame,
    domain: Domain,
    budget: Interval,
}

impl FaceDomain {
    fn new(context: &BandContext<'_>, face: &Face, budget: f64) -> Option<Self> {
        let wires = face.boundaries();
        let wire = wires.first()?;
        let construct = |frame, positive| {
            Some(Self {
                frame,
                domain: Domain::new(context, frame, wire, positive)?,
                budget: chart_budget(frame, budget)?,
            })
        };
        if let Some(frame) = Frame::new(
            context,
            Ring {
                axis: 0,
                fixed: 0.0,
                positive: false,
            },
            1.0,
        ) && let Some(domain) = construct(frame, context.outward)
        {
            return Some(domain);
        }
        match context.carrier {
            SurfaceCarrier::Plane { origin, normal } => construct(
                projection_frame(context, origin, normal, false),
                context.outward,
            ),
            SurfaceCarrier::LinearRevolution {
                origin,
                axis,
                slope,
                ..
            } if slope != 0.0 => construct(
                projection_frame(context, origin, axis, false),
                context.outward == (slope < 0.0),
            ),
            SurfaceCarrier::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                ..
            } if radial_center.abs() <= crate::KERNEL_EPSILON => {
                let sum = wire.iter().fold(Vector3::new(0.0, 0.0, 0.0), |sum, edge| {
                    sum + (edge.front().point() - axial_center)
                });
                let direction = (sum.magnitude2() > 0.0 && sum.magnitude2().is_finite())
                    .then(|| sum.normalize());
                direction
                    .into_iter()
                    .chain([axis, -axis])
                    .find_map(|normal| {
                        construct(
                            projection_frame(context, axial_center, normal, true),
                            context.outward,
                        )
                    })
            }
            _ => None,
        }
    }

    fn excludes(&self, point: Point3) -> bool {
        if matches!(self.frame.chart, Chart::Orthographic { hemisphere: true }) {
            let [_, _, z, _] = self
                .frame
                .local([point.x, point.y, point.z, 1.0].map(Interval::point));
            if z.add(self.budget).hi < 0.0 {
                return true;
            }
        }
        source_ball(self.frame, point, self.budget)
            .is_some_and(|query| self.domain.excludes(self.frame, query))
    }
}

pub(in super::super::super::super) fn missing_source_vertex(
    context: &BandContext<'_>,
    faces: &[Face],
    budget: f64,
) -> Option<usize> {
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
    // Existing topological endpoints are explicit positive witnesses.
    // Avoid constructing interval charts for ordinary unchanged facets.
    let key = |p: Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let mut boundary_points = BTreeSet::new();
    for face in faces {
        for wire in face.boundaries() {
            boundary_points.extend(wire.iter().map(|edge| key(edge.front().point())));
        }
    }
    vertices.retain(|&id| {
        context
            .evidence
            .positions
            .get(id)
            .is_some_and(|&point| !boundary_points.contains(&key(point)))
    });
    if vertices.is_empty() {
        return None;
    }
    // Each chart face may choose a different coordinate system. All must
    // exclude the source ball before rejecting the region. An unresolved
    // face never counts as an exclusion.
    let domains = faces
        .iter()
        .map(|face| FaceDomain::new(context, face, budget))
        .collect::<Option<Vec<_>>>()?;
    if domains.is_empty() {
        return None;
    }
    vertices.into_iter().find(|&id| {
        context
            .evidence
            .positions
            .get(id)
            .is_some_and(|&point| domains.iter().all(|domain| domain.excludes(point)))
    })
}

#[cfg(test)]
mod tests {
    use super::super::Controls;
    use super::*;

    fn frame() -> Frame {
        Frame {
            chart: Chart::Revolution,
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
