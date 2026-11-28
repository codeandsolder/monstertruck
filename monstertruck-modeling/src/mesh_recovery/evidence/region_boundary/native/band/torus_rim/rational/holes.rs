//! Certify bounded, negatively oriented contractible holes in an annular chart.
use super::{BandContext, Edge, Frame, Interval, Rejection, Ring, lift_cycle, simplicity};
use std::f64::consts::TAU;

fn signed_area(pieces: &[super::Piece], origin: f64) -> Interval {
    // Green's theorem: integral x dy on a closed loop. Each piece's
    // parameter is normalized to [0,1], so its rate bounds integrate
    // directly. The translation by origin improves cancellation.
    pieces.iter().fold(Interval::point(0.0), |sum, p| {
        sum.add(p.bounds[0].sub(Interval::point(origin)).mul(p.rates[1]))
    })
}

fn oriented_hole(
    frame: Frame,
    mut pieces: Vec<super::Piece>,
    origin: f64,
    expected: f64,
) -> Option<()> {
    for _ in 0..=8 {
        if signed_area(&pieces, origin).scale(expected).lo > 0.0 {
            return Some(());
        }
        if pieces.len() > super::MAX_PIECES / 2 {
            return None;
        }
        pieces = pieces
            .iter()
            .map(|p| p.split(frame))
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect();
    }
    None
}

pub(in super::super::super) fn hole_bounds(
    context: &BandContext<'_>,
    rim: &[Edge],
    axis: usize,
    lower: f64,
    span: f64,
) -> Result<[[f64; 2]; 2], Rejection> {
    if rim.is_empty()
        || rim
            .iter()
            .zip(rim.iter().cycle().skip(1))
            .any(|(a, b)| a.back() != b.front())
    {
        return Err(Rejection::Invalid);
    }
    let frame = Frame::new(
        context,
        Ring {
            axis,
            fixed: lower,
            positive: false,
        },
        1.0,
    )
    .ok_or(Rejection::Invalid)?;
    let uv = context
        .torus_chart(rim.first().ok_or(Rejection::Invalid)?.front().point())
        .ok_or(Rejection::Invalid)?;
    let start = [uv[axis], lower + (uv[1 - axis] - lower).rem_euclid(TAU)].map(Interval::point);
    let (pieces, end) = lift_cycle(frame, rim, start)?;
    let resolution = frame.resolution;
    if (0..2).any(|i| end[i].sub(start[i]).norm() > resolution) {
        return Err(Rejection::Winding);
    }
    let bounds = pieces
        .iter()
        .map(|p| p.bounds)
        .reduce(|a, b| [a[0].union(b[0]), a[1].union(b[1])])
        .ok_or(Rejection::Invalid)?;
    if bounds[0].width() >= TAU - resolution
        || bounds[1].lo <= lower + resolution
        || bounds[1].hi >= lower + span - resolution
    {
        return Err(Rejection::Invalid);
    }
    simplicity::cycle(frame, &pieces, false).map_err(Rejection::Pair)?;
    let expected = if context.outward == (axis == 0) {
        -1.0
    } else {
        1.0
    };
    oriented_hole(frame, pieces, bounds[0].mid(), expected).ok_or(Rejection::Invalid)?;
    Ok(bounds.map(|x| [x.lo, x.hi]))
}
