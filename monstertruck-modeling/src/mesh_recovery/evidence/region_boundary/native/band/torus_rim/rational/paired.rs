//! Two variable rims on a regular revolution annulus. Prove the actual
//! rational curves are simple, oppositely winding and mutually disjoint in
//! a common periodic strip, then cut both native cycles at the same two
//! unique transverse coordinates. No reference-ring approximation is used.

use super::super::super::{
    BandPlan, NativeChartObstacle, NativeTrimError, chart_point, commit_rings, copy_rings, cut_ring,
};
use super::super::{OrientedRim, chart_obstacles, regular_radial_band, signed_delta};
use super::{
    BandContext, Frame, Interval, Piece, Rejection, Ring, lift_cycle, near, simplicity, unique_cut,
};
use crate::mesh_recovery::{OrientedTrim, SurfaceCarrier};
use crate::{Edge, InnerSpace, ParametricCurve, Point3, Vertex};
use std::f64::consts::TAU;

struct LiftedRim {
    rim: OrientedRim,
    pieces: Vec<Piece>,
    height: Interval,
}

#[derive(Clone, Copy)]
struct Cut {
    angle: f64,
    points: [Point3; 2],
    heights: [f64; 2],
}

impl LiftedRim {
    fn new(context: &BandContext<'_>, frame: Frame, rim: OrientedRim) -> Result<Self, Rejection> {
        let uv = chart_point(
            context.carrier,
            rim.edges.first().ok_or(Rejection::Invalid)?.front().point(),
            context.options.angular_resolution,
        )
        .ok_or(Rejection::Invalid)?;
        let start = [uv[frame.varying], uv[1 - frame.varying]].map(Interval::point);
        let (pieces, end) = lift_cycle(frame, &rim.edges, start)?;
        if !frame.junction(end, start, TAU) {
            return Err(Rejection::Winding);
        }
        let height = pieces
            .iter()
            .map(|p| p.bounds[1])
            .reduce(Interval::union)
            .ok_or(Rejection::Invalid)?;
        if frame.periodic_height() && height.width() >= TAU - frame.height_resolution {
            return Err(Rejection::Invalid);
        }
        Ok(Self {
            rim,
            pieces,
            height,
        })
    }

    fn shift_height(&mut self, shift: f64) {
        let shift = Interval::point(shift);
        self.height = self.height.add(shift);
        for piece in &mut self.pieces {
            piece.start[1] = piece.start[1].add(shift);
            piece.end[1] = piece.end[1].add(shift);
            piece.bounds[1] = piece.bounds[1].add(shift);
        }
    }

    fn point_phase(
        &self,
        context: &BandContext<'_>,
        frame: Frame,
        index: usize,
        t: f64,
    ) -> Option<(Point3, [f64; 2])> {
        let piece = &self.pieces[index];
        let point = self.rim.edges[piece.edge].oriented_curve().evaluate(t);
        let uv = chart_point(context.carrier, point, context.options.angular_resolution)?;
        let phase = [uv[frame.varying], uv[1 - frame.varying]];
        Some((
            point,
            std::array::from_fn(|i| {
                frame
                    .unwrap(i, Interval::point(phase[i]), piece.bounds[i].mid())
                    .mid()
            }),
        ))
    }

    fn interior(&self, context: &BandContext<'_>, piece: &Piece, point: Point3) -> bool {
        let edge = &self.rim.edges[piece.edge];
        let clearance = 4.0 * context.options.intersection_tolerance;
        (point - edge.front().point()).magnitude() > clearance
            && (point - edge.back().point()).magnitude() > clearance
    }

    fn candidate(
        &self,
        context: &BandContext<'_>,
        frame: Frame,
        index: usize,
    ) -> Option<(Point3, [f64; 2])> {
        let piece = &self.pieces[index];
        if piece.rates[0].lo <= 0.0 {
            return None;
        }
        let (point, phase) = self.point_phase(
            context,
            frame,
            index,
            piece.controls.parameter[0].midpoint(piece.controls.parameter[1]),
        )?;
        (self.interior(context, piece, point)
            && phase[0] > piece.start[0].hi + frame.resolution
            && phase[0] < piece.end[0].lo - frame.resolution
            && unique_cut(&self.pieces, index, phase[0], frame.resolution))
        .then_some((point, phase))
    }

    // The target phase has exactly one intersection with the entire cycle,
    // inside a continuously certified strictly advancing rational span.
    // Bisection only locates that already proved root on its actual curve.
    fn cut_at(&self, context: &BandContext<'_>, frame: Frame, angle: f64) -> Option<(Point3, f64)> {
        let mut matches = self.pieces.iter().enumerate().filter_map(|(i, p)| {
            let target = near(Interval::point(angle), p.bounds[0].mid()).mid();
            (p.rates[0].lo > 0.0
                && target > p.start[0].hi + frame.resolution
                && target < p.end[0].lo - frame.resolution)
                .then_some((i, target))
        });
        let (index, target) = matches.next()?;
        if matches.next().is_some() || !unique_cut(&self.pieces, index, target, frame.resolution) {
            return None;
        }
        let piece = &self.pieces[index];
        let [mut lo, mut hi] = piece.controls.parameter;
        for _ in 0..64 {
            let mid = lo.midpoint(hi);
            if mid <= lo || mid >= hi {
                break;
            }
            let (_, phase) = self.point_phase(context, frame, index, mid)?;
            if phase[0] < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let (point, phase) = self.point_phase(context, frame, index, lo.midpoint(hi))?;
        (self.interior(context, piece, point) && (phase[0] - target).abs() <= frame.resolution)
            .then_some((point, phase[1]))
    }
}

fn find_cuts(context: &BandContext<'_>, frame: Frame, rims: [&LiftedRim; 2]) -> Option<[Cut; 2]> {
    let mut candidates = Vec::new();
    for i in 0..rims[0].pieces.len() {
        let Some((a, phase)) = rims[0].candidate(context, frame, i) else {
            continue;
        };
        let Some((b, height)) = rims[1].cut_at(context, frame, phase[0]) else {
            continue;
        };
        candidates.push(Cut {
            angle: phase[0],
            points: [a, b],
            heights: [phase[1], height],
        });
    }
    let first = *candidates.first()?;
    let distance = |p: &Cut| signed_delta(first.angle, p.angle).abs();
    let second = *candidates
        .iter()
        .max_by(|a, b| distance(a).total_cmp(&distance(b)))?;
    (distance(&second) > frame.resolution).then_some([first, second])
}

struct PairProof {
    axis: usize,
    cuts: [Cut; 2],
}

fn prove_pair(
    context: &BandContext<'_>,
    rims: &[OrientedRim; 2],
    axis: usize,
    lower: usize,
    obstacles: &mut Vec<NativeChartObstacle>,
) -> Option<PairProof> {
    let frame = Frame::new(
        context,
        Ring {
            axis,
            fixed: 0.0,
            positive: false,
        },
        1.0,
    )?;
    // For an outward chart, positive azimuth winding has material on its
    // upper meridian side; positive meridional winding has material on
    // its lower azimuth side. Reverse the negative cycle for the proof.
    let lower_positive = context.outward == (axis == 0);
    let mut lift = |rim: OrientedRim| {
        // Keep the source mapping even when constructing the lift fails.
        let sources = rim.sources.clone();
        match LiftedRim::new(context, frame, rim) {
            Ok(lifted) => Some(lifted),
            Err(error) => {
                obstacles.extend(chart_obstacles(error, &sources));
                None
            }
        }
    };
    let a = lift(rims[lower].reversed(!lower_positive))?;
    let mut b = lift(rims[1 - lower].reversed(lower_positive))?;
    for rim in [&a, &b] {
        if let Err(pairs) = simplicity::cycle(frame, &rim.pieces, true) {
            obstacles.extend(chart_obstacles(Rejection::Pair(pairs), &rim.rim.sources));
            return None;
        }
    }
    let mut cuts = find_cuts(context, frame, [&a, &b]).or_else(|| {
        obstacles.push(NativeChartObstacle::CutSelection);
        None
    })?;
    // Fix the transverse lift using the observed oriented cut, not the
    // relative positions of independently chosen angular branch origins.
    if frame.periodic_height() {
        let [lo, hi] = cuts[0].heights;
        let shift = ((lo - hi) / TAU).floor().mul_add(TAU, TAU);
        b.shift_height(shift);
        for cut in &mut cuts {
            cut.heights[1] += shift;
        }
    }
    let extent = a.height.union(b.height);
    let regular = match context.carrier {
        SurfaceCarrier::CircularRevolution {
            radial_center,
            profile_radius,
            ..
        } => {
            extent.width() < TAU - frame.height_resolution
                && regular_radial_band(
                    radial_center,
                    profile_radius,
                    axis,
                    extent.lo,
                    extent.width(),
                    context.options.intersection_tolerance,
                )
        }
        SurfaceCarrier::LinearRevolution {
            radius_at_origin,
            slope,
            ..
        } => {
            slope
                .mul_add(extent.lo, radius_at_origin)
                .min(slope.mul_add(extent.hi, radius_at_origin))
                > context.options.intersection_tolerance
        }
        SurfaceCarrier::Plane { .. } => false,
    };
    if !regular
        || cuts
            .iter()
            .any(|c| c.heights[1] - c.heights[0] <= frame.height_resolution)
        || !context.covers_band(1 - axis, extent.lo, extent.width())
        || !context.source_chart_positive()
    {
        obstacles.push(NativeChartObstacle::Coverage);
        return None;
    }
    if let Err(pairs) = simplicity::disjoint_cycles(frame, &a.pieces, &b.pieces) {
        obstacles.extend(pairs.into_iter().filter_map(|[i, j]| {
            Some(NativeChartObstacle::RimPair([
                *a.rim.sources.get(i)?,
                *b.rim.sources.get(j)?,
            ]))
        }));
        return None;
    }
    // Put both cut endpoints back in original wire order. Both seams run
    // from wire zero to wire one, as required by the shared BandPlan.
    if lower == 1 {
        for cut in &mut cuts {
            cut.points.swap(0, 1);
            cut.heights.swap(0, 1);
        }
    }
    Some(PairProof { axis, cuts })
}

impl BandContext<'_> {
    pub(in super::super::super) fn paired_band(
        &self,
        cycles: [&[OrientedTrim]; 2],
        trims: &mut [Result<Vec<Edge>, NativeTrimError>],
        obstacles: &mut Vec<NativeChartObstacle>,
    ) -> Option<BandPlan> {
        let mut axes = match self.carrier {
            SurfaceCarrier::LinearRevolution { .. } => 0..1,
            SurfaceCarrier::CircularRevolution { radial_center, .. }
                if radial_center > crate::KERNEL_EPSILON =>
            {
                0..2
            }
            _ => return None,
        };
        let mut changed = copy_rings(cycles[0].iter().chain(cycles[1]), trims)?;
        let rims = [
            OrientedRim::collect(cycles[0], |id| changed.get(&id).map(Vec::as_slice))?,
            OrientedRim::collect(cycles[1], |id| changed.get(&id).map(Vec::as_slice))?,
        ];
        let proof = axes.find_map(|axis| {
            (0..2).find_map(|lower| prove_pair(self, &rims, axis, lower, obstacles))
        })?;
        let plan = (|| {
            let mut vertices = [Vec::<Vertex>::new(), Vec::new()];
            let mut seams = Vec::new();
            for cut in proof.cuts {
                for i in 0..2 {
                    vertices[i].push(cut_ring(
                        cycles[i],
                        cut.points[i],
                        &mut changed,
                        self.options.intersection_tolerance,
                    )?);
                }
                seams.push(self.seam(
                    proof.axis,
                    cut.angle,
                    cut.heights[0],
                    cut.heights[1] - cut.heights[0],
                    vertices[0].last()?,
                    vertices[1].last()?,
                )?);
            }
            let [Some(a), Some(b)] = vertices.map(|v| v.try_into().ok()) else {
                return None;
            };
            Some(BandPlan::unperforated([a, b], seams.try_into().ok()?))
        })();
        if plan.is_some() {
            commit_rings(changed, trims);
        } else {
            obstacles.push(NativeChartObstacle::SeamConstruction);
        }
        plan
    }
}
