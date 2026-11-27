//! Continuous certificate for the actual rational native rim, including
//! polynomial approximations of regular analytic intersection branches.
//! No tensor-circle assumption is made for a general section spline.
mod controls;
mod holes;
pub(in super::super) use holes::hole_bounds;
mod interval;
mod simplicity;
#[cfg(test)]
mod tests;

use super::{BandContext, Ring, signed_delta};
use crate::mesh_recovery::SurfaceCarrier;
use crate::{Edge, InnerSpace, ParametricCurve, Point3, Vector3};
use controls::{Control, Controls};
use interval::{Interval, angle};
use monstertruck_geometry::prelude::TryIntoHomogeneousBsplineCurve;
use std::f64::consts::{PI, TAU};

const MAX_PIECES: usize = 8192;
const MAX_DEPTH: usize = 16;

#[derive(Clone, Copy)]
struct Frame {
    center: Point3,
    basis: [Vector3; 3],
    major: f64,
    varying: usize,
    direction: f64,
    side: f64,
    resolution: f64,
}

impl Frame {
    fn new(context: &BandContext<'_>, ring: Ring, side: f64) -> Option<Self> {
        let SurfaceCarrier::CircularRevolution {
            axial_center: center,
            axis,
            radial_center: major,
            ..
        } = context.carrier
        else {
            return None;
        };
        let (u, v) = super::super::normal_basis(axis);
        Some(Self {
            center,
            basis: [u, v, axis],
            major,
            varying: ring.axis,
            direction: if ring.positive { -1.0 } else { 1.0 },
            side,
            resolution: context.options.angular_resolution,
        })
    }
    fn local(self, value: Control) -> Control {
        let centered: [Interval; 3] =
            std::array::from_fn(|i| value[i].sub(value[3].scale(self.center[i])));
        let mut result = value;
        for (i, basis) in self.basis.iter().enumerate() {
            result[i] = (0..3)
                .map(|j| centered[j].scale(basis[j]))
                .fold(Interval::point(0.0), Interval::add);
        }
        result
    }
    fn orient(self, coordinates: [Interval; 2]) -> [Interval; 2] {
        [
            coordinates[self.varying].scale(self.direction),
            coordinates[1 - self.varying].scale(self.side),
        ]
    }
    fn endpoint(self, control: Control) -> Option<[Interval; 2]> {
        let [x, y, z, w] = control;
        let [x, y, z] = [x.div(w)?, y.div(w)?, z.div(w)?];
        let rho = x.square().add(y.square()).sqrt()?;
        Some(self.orient([
            angle(x, y)?,
            angle(rho.sub(Interval::point(self.major)), z)?,
        ]))
    }
    fn rates(self, controls: &[Control]) -> Option<[Interval; 2]> {
        let hull = |i: usize| controls.iter().map(|p| p[i]).reduce(Interval::union);
        let weight = hull(3)?;
        if weight.lo <= 0.0 {
            return None;
        }
        let degree = crate::usize_to_f64(controls.len() - 1);
        let derivative = |i: usize| {
            controls
                .windows(2)
                .map(|p| p[1][i].sub(p[0][i]).scale(degree))
                .reduce(Interval::union)
        };
        let weight_rate = derivative(3)?;
        let cartesian = |i: usize| -> Option<Interval> {
            let mut bounds = controls.first()?.get(i)?.div(controls.first()?[3])?;
            for point in &controls[1..] {
                bounds = bounds.union(point[i].div(point[3])?);
            }
            Some(bounds)
        };
        let [x, y, z] = [cartesian(0)?, cartesian(1)?, cartesian(2)?];
        let rate = |i: usize, p: Interval| derivative(i)?.sub(p.mul(weight_rate)).div(weight);
        let [dx, dy, dz] = [rate(0, x)?, rate(1, y)?, rate(2, z)?];
        let radius_squared = x.square().add(y.square());
        let radius = radius_squared.sqrt()?;
        if radius.lo <= self.resolution {
            return None;
        }
        let radial_rate = x.mul(dx).add(y.mul(dy)).div(radius)?;
        let profile_radial = radius.sub(Interval::point(self.major));
        let meridian2 = profile_radial.square().add(z.square());
        if meridian2.lo <= self.resolution * self.resolution {
            return None;
        }
        Some(
            self.orient([
                x.mul(dy).sub(y.mul(dx)).div(radius_squared)?,
                profile_radial
                    .mul(dz)
                    .sub(z.mul(radial_rate))
                    .div(meridian2)?,
            ]),
        )
    }
}

// Unwrap only within a proved short local chart. The winding accumulator
// uses absolute angular lifts rather than summing approximate increments.
fn near(angle: Interval, reference: f64) -> Interval {
    let shift = ((reference - angle.mid()) / TAU).round() * TAU;
    angle.add(Interval::point(shift))
}

struct Piece {
    controls: Controls,
    start: [Interval; 2],
    end: [Interval; 2],
    bounds: [Interval; 2],
    rates: [Interval; 2],
    edge: usize,
}

impl Piece {
    fn new(
        frame: Frame,
        controls: Controls,
        reference: [Interval; 2],
        edge: usize,
    ) -> Option<Self> {
        let rates = frame.rates(&controls.values)?;
        if rates
            .iter()
            .any(|r| !r.finite() || r.magnitude() >= PI * 0.5)
            || !rates.iter().any(|r| r.directed())
        {
            return None;
        }
        let raw_start = frame.endpoint(*controls.values.first()?)?;
        let start = std::array::from_fn(|i| near(raw_start[i], reference[i].mid()));
        if (0..2).any(|i| start[i].sub(reference[i]).magnitude() > frame.resolution) {
            return None;
        }
        let raw_end = frame.endpoint(*controls.values.last()?)?;
        let end = std::array::from_fn(|i| near(raw_end[i], start[i].mid()));
        let extent = |i: usize| {
            let step = Interval::point(0.0).union(rates[i]);
            start[i].add(step).intersection(end[i].sub(step))
        };
        let bounds = [extent(0)?, extent(1)?];
        Some(Self {
            controls,
            start,
            end,
            bounds,
            rates,
            edge,
        })
    }
    fn shifted_bounds(&self, shift: f64) -> [Interval; 2] {
        [self.bounds[0].add(Interval::point(shift)), self.bounds[1]]
    }
    fn split(&self, frame: Frame) -> Option<[Self; 2]> {
        let [a, b] = self.controls.split();
        let first = Self::new(frame, a, self.start, self.edge)?;
        let second = Self::new(frame, b, first.end, self.edge)?;
        Some([first, second])
    }
}

fn append(
    frame: Frame,
    controls: &Controls,
    reference: [Interval; 2],
    edge: usize,
    depth: usize,
    pieces: &mut Vec<Piece>,
) -> Option<[Interval; 2]> {
    if pieces.len() >= MAX_PIECES {
        return None;
    }
    // Retain the interval polynomial for fallback subdivision.
    let copy = Controls {
        values: controls.values.clone(),
        parameter: controls.parameter,
    };
    if let Some(piece) = Piece::new(frame, copy, reference, edge) {
        let end = piece.end;
        pieces.push(piece);
        return Some(end);
    }
    if depth == 0 {
        return None;
    }
    let [a, b] = controls.split();
    let mid = append(frame, &a, reference, edge, depth - 1, pieces)?;
    append(frame, &b, mid, edge, depth - 1, pieces)
}

pub(super) struct RimProof {
    pub(super) extent: f64,
    pub(super) cuts: [([f64; 2], Point3); 2],
}

fn unique_cut(pieces: &[Piece], index: usize, angle: f64, resolution: f64) -> bool {
    pieces.iter().enumerate().all(|(i, p)| {
        if i == index {
            return true;
        }
        let phase = near(Interval::point(angle), p.bounds[0].mid());
        p.bounds[0].sub(phase).lo > resolution || p.bounds[0].sub(phase).hi < -resolution
    })
}

#[derive(Debug)]
pub(in super::super) enum Rejection {
    Span(usize),
    Winding,
    Ring(usize),
    Pair(Vec<[usize; 2]>),
    Cuts,
    Invalid,
}

fn lift_cycle(
    frame: Frame,
    rim: &[Edge],
    first: [Interval; 2],
) -> Result<(Vec<Piece>, [Interval; 2]), Rejection> {
    let mut reference = first;
    let mut pieces = Vec::new();
    for (edge_id, edge) in rim.iter().enumerate() {
        let curve = edge
            .oriented_curve()
            .try_into_homogeneous_bspline_curve()
            .ok_or(Rejection::Span(edge_id))?;
        for mut span in controls::spans(&curve).ok_or(Rejection::Span(edge_id))? {
            for control in &mut span.values {
                *control = frame.local(*control);
            }
            reference = append(frame, &span, reference, edge_id, MAX_DEPTH, &mut pieces)
                .ok_or(Rejection::Span(edge_id))?;
        }
    }
    Ok((pieces, reference))
}

pub(super) fn prove(
    context: &BandContext<'_>,
    rim: &[Edge],
    ring: Ring,
    side: f64,
) -> Result<RimProof, Rejection> {
    let frame = Frame::new(context, ring, side).ok_or(Rejection::Invalid)?;
    let uv = context
        .torus_chart(rim.first().ok_or(Rejection::Invalid)?.front().point())
        .ok_or(Rejection::Invalid)?;
    let first_height = (side * (uv[1 - ring.axis] - ring.fixed)).rem_euclid(TAU);
    let first = [
        Interval::point(frame.direction * uv[ring.axis]),
        Interval::point(side.mul_add(ring.fixed, first_height)),
    ];
    let (pieces, reference) = lift_cycle(frame, rim, first)?;
    let resolution = frame.resolution;
    if reference[0]
        .sub(first[0])
        .sub(Interval::point(TAU))
        .magnitude()
        > resolution
        || reference[1].sub(first[1]).magnitude() > resolution
    {
        return Err(Rejection::Winding);
    }
    let height_origin = side * ring.fixed;
    let mut extent = 0.0_f64;
    for piece in &pieces {
        let height = piece.bounds[1].sub(Interval::point(height_origin));
        if height.lo <= resolution || height.hi >= TAU - resolution {
            return Err(Rejection::Ring(piece.edge));
        }
        extent = extent.max(height.hi);
    }
    simplicity::cycle(frame, &pieces, true).map_err(Rejection::Pair)?;
    let mut candidates = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        if piece.rates[0].lo <= 0.0 {
            continue;
        }
        let curve = rim[piece.edge].oriented_curve();
        let point =
            curve.evaluate(piece.controls.parameter[0].midpoint(piece.controls.parameter[1]));
        let clearance = 4.0 * context.options.intersection_tolerance;
        if (point - rim[piece.edge].front().point()).magnitude() <= clearance
            || (point - rim[piece.edge].back().point()).magnitude() <= clearance
        {
            continue;
        }
        let uv = context.torus_chart(point).ok_or(Rejection::Cuts)?;
        let phase = near(
            Interval::point(frame.direction * uv[ring.axis]),
            piece.bounds[0].mid(),
        )
        .mid();
        if phase <= piece.start[0].hi + resolution
            || phase >= piece.end[0].lo - resolution
            || !unique_cut(&pieces, i, phase, resolution)
        {
            continue;
        }
        candidates.push((uv, point));
    }
    let first = *candidates.first().ok_or(Rejection::Cuts)?;
    let distance = |p: &([f64; 2], Point3)| signed_delta(first.0[ring.axis], p.0[ring.axis]).abs();
    let second = *candidates
        .iter()
        .max_by(|a, b| distance(a).total_cmp(&distance(b)))
        .ok_or(Rejection::Cuts)?;
    (distance(&second) > resolution)
        .then_some(RimProof {
            extent,
            cuts: [first, second],
        })
        .ok_or(Rejection::Cuts)
}
