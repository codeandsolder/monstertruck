//! Periodic simple-cycle proof on actual rational chart enclosures.
use super::{Frame, Piece, interval::Interval};
use std::f64::consts::TAU;

#[derive(Clone, Copy)]
enum Junction {
    EndStart,
    StartEnd,
}

fn separate(a: [Interval; 2], b: [Interval; 2], margin: f64) -> bool {
    (0..2).any(|i| a[i].hi + margin < b[i].lo || b[i].hi + margin < a[i].lo)
}

fn tangent_separates(a: &Piece, b: &Piece, shift: f64, junction: Junction, margin: f64) -> bool {
    let (left, right, offset) = match junction {
        Junction::EndStart => (a, b, shift),
        Junction::StartEnd => (b, a, -shift),
    };
    if left.end[0]
        .sub(right.start[0].add(Interval::point(offset)))
        .magnitude()
        > margin
        || left.end[1].sub(right.start[1]).magnitude() > margin
    {
        return false;
    }
    let normalized = |p: &Piece| {
        let [x, y] = p.rates.map(Interval::mid);
        let len = x.hypot(y);
        (len.is_finite() && len > 0.0).then_some([x / len, y / len])
    };
    let (Some(a), Some(b)) = (normalized(left), normalized(right)) else {
        return false;
    };
    let d = [a[0] + b[0], a[1] + b[1]];
    [left, right]
        .iter()
        .all(|p| p.rates[0].scale(d[0]).add(p.rates[1].scale(d[1])).lo > 0.0)
}

fn disjoint(
    frame: Frame,
    a: &Piece,
    b: &Piece,
    shift: f64,
    junction: Option<Junction>,
    depth: usize,
    remaining: &mut usize,
) -> bool {
    if *remaining == 0 {
        return false;
    }
    *remaining -= 1;
    let bounds = b.shifted_bounds(shift);
    let margin = frame.resolution;
    if separate(a.bounds, bounds, margin)
        || junction.is_some_and(|j| tangent_separates(a, b, shift, j, margin))
    {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let extent = |p: &Piece| p.bounds.iter().map(|b| b.width()).fold(0.0_f64, f64::max);
    if extent(a) >= extent(b) {
        let Some([first, last]) = a.split(frame) else {
            return false;
        };
        disjoint(
            frame,
            &first,
            b,
            shift,
            junction.filter(|j| matches!(j, Junction::StartEnd)),
            depth - 1,
            remaining,
        ) && disjoint(
            frame,
            &last,
            b,
            shift,
            junction.filter(|j| matches!(j, Junction::EndStart)),
            depth - 1,
            remaining,
        )
    } else {
        let Some([first, last]) = b.split(frame) else {
            return false;
        };
        disjoint(
            frame,
            a,
            &first,
            shift,
            junction.filter(|j| matches!(j, Junction::EndStart)),
            depth - 1,
            remaining,
        ) && disjoint(
            frame,
            a,
            &last,
            shift,
            junction.filter(|j| matches!(j, Junction::StartEnd)),
            depth - 1,
            remaining,
        )
    }
}

// The periodic chart boxes are shorter than a half-turn. Place every
// box's midpoint in one period and include its two nearest images.
// Compare only overlapping intervals, with the smaller-index piece in
// its base image, so every potentially overlapping pair appears once.
struct Image {
    piece: usize,
    shift: f64,
    base: bool,
    interval: Interval,
}

fn sweep_images(pieces: &[Piece]) -> Vec<Image> {
    let mut images = Vec::with_capacity(3 * pieces.len());
    for (piece, p) in pieces.iter().enumerate() {
        let period = (p.bounds[0].mid() / TAU).floor();
        for image in [-1.0, 0.0, 1.0] {
            let shift = (image - period) * TAU;
            images.push(Image {
                piece,
                shift,
                base: image == 0.0,
                interval: p.bounds[0].add(Interval::point(shift)),
            });
        }
    }
    images.sort_unstable_by(|a, b| a.interval.lo.total_cmp(&b.interval.lo));
    images
}

fn junction(i: usize, j: usize, count: usize, shift: f64) -> Option<Junction> {
    if j == i + 1 && shift.abs() < TAU * 0.5 {
        Some(Junction::EndStart)
    } else if i == 0 && j + 1 == count && (shift + TAU).abs() < TAU * 0.5 {
        Some(Junction::StartEnd)
    } else {
        None
    }
}

pub(super) fn cycle(frame: Frame, pieces: &[Piece]) -> Result<(), Vec<[usize; 2]>> {
    if pieces.len() < 2 {
        return Err(Vec::new());
    }
    // The caller has already certified a continuous winding-one lift.
    // Strict monotonicity in its periodic coordinate alone proves a
    // simple cycle, without requiring quadratic pair comparisons.
    if pieces.iter().all(|p| p.rates[0].lo > 0.0) {
        return Ok(());
    }
    let images = sweep_images(pieces);
    let mut remaining = 1 << 20;
    let mut rejected = std::collections::BTreeSet::new();
    for (k, image) in images.iter().enumerate() {
        for other in &images[k + 1..] {
            if other.interval.lo > image.interval.hi + frame.resolution {
                break;
            }
            let (left, right) = if image.piece < other.piece {
                (image, other)
            } else {
                (other, image)
            };
            if left.piece == right.piece || !left.base {
                continue;
            }
            let [a, b] = [&pieces[left.piece], &pieces[right.piece]];
            let shift = right.shift - left.shift;
            let adjacent = junction(left.piece, right.piece, pieces.len(), shift);
            let pair = [a.edge.min(b.edge), a.edge.max(b.edge)];
            if rejected.contains(&pair) {
                continue;
            }
            if !disjoint(frame, a, b, shift, adjacent, 24, &mut remaining) {
                rejected.insert(pair);
                if remaining == 0 || rejected.len() >= 64 {
                    return Err(rejected.into_iter().collect());
                }
            }
        }
    }
    if rejected.is_empty() {
        Ok(())
    } else {
        Err(rejected.into_iter().collect())
    }
}
