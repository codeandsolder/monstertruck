//! Shared best-first Bernstein span subdivision for global deviation brackets.
//!
//! Numerical conditioning allowances belong to each carrier-specific caller.

use crate::analytic_bezier_partition::ProofBezierSpan;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

const MAX_SUBDIVISIONS: usize = 32_768;
const MAX_SUBDIVISION_DEPTH: usize = 48;

#[derive(Clone, Copy)]
pub struct MaximumBracket {
    pub lower: f64,
    pub upper: f64,
    pub subdivisions: usize,
    pub max_depth: usize,
}

struct Span {
    controls: Vec<[f64; 4]>,
    upper: f64,
    depth: usize,
    sequence: usize,
}

impl PartialEq for Span {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.upper.total_cmp(&other.upper).is_eq()
    }
}

impl Eq for Span {}

impl PartialOrd for Span {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}

impl Ord for Span {
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper
            .total_cmp(&other.upper)
            .then(self.sequence.cmp(&other.sequence))
    }
}

type Controls = Vec<[f64; 4]>;
type Split = (Controls, Controls);

/// Homogeneous de Casteljau subdivision at the exact parameter midpoint.
fn subdivide(controls: &[[f64; 4]]) -> Option<Split> {
    let mut row = controls.to_vec();
    let mut left = Vec::with_capacity(row.len());
    let mut right = Vec::with_capacity(row.len());
    while let (Some(first), Some(last)) = (row.first(), row.last()) {
        left.push(*first);
        right.push(*last);
        if row.len() == 1 {
            break;
        }
        row = row
            .windows(2)
            .map(|pair| std::array::from_fn(|i| f64::midpoint(pair[0][i], pair[1][i])))
            .collect();
        if row.iter().flatten().any(|value| !value.is_finite()) {
            return None;
        }
    }
    right.reverse();
    Some((left, right))
}

/// Refine global maximum bounds for each original positive knot span.
///
/// The evaluator supplies (attained endpoint witness, global Bernstein upper
/// bound) for the *whole* supplied Bézier subspan. The algorithm refines the
/// largest bound first, retaining a valid but potentially wider bracket if the
/// work/depth limit is reached. This uses f64 control arithmetic: callers must
/// separately pad the bracket for floating-point conditioning.
pub fn refine_maximum(
    spans: Vec<ProofBezierSpan>,
    target_gap: f64,
    evaluate: impl Fn(&[[f64; 4]]) -> Option<(f64, f64)>,
) -> Option<MaximumBracket> {
    if !target_gap.is_finite() || target_gap <= 0.0 || spans.is_empty() {
        return None;
    }
    let mut nodes = BinaryHeap::with_capacity(spans.len());
    let mut lower = 0.0_f64;
    for (sequence, span) in spans.into_iter().enumerate() {
        let (witness, upper) = evaluate(&span.controls)?;
        if !witness.is_finite() || !upper.is_finite() || witness < 0.0 || upper < witness {
            return None;
        }
        lower = lower.max(witness);
        nodes.push(Span {
            controls: span.controls,
            upper,
            depth: 0,
            sequence,
        });
    }
    let mut next_sequence = nodes.len();
    let mut max_depth = 0;
    let mut subdivisions = 0;
    loop {
        let largest = nodes.peek()?;
        if largest.upper - lower <= target_gap
            || largest.depth >= MAX_SUBDIVISION_DEPTH
            || subdivisions >= MAX_SUBDIVISIONS
        {
            break;
        }
        let previous = nodes.pop()?;
        let depth = previous.depth + 1;
        let (left, right) = subdivide(&previous.controls)?;
        let children: [Controls; 2] = (left, right).into();
        for controls in children {
            let (witness, upper) = evaluate(&controls)?;
            if !witness.is_finite() || !upper.is_finite() || witness < 0.0 || upper < witness {
                return None;
            }
            lower = lower.max(witness);
            nodes.push(Span {
                controls,
                upper,
                depth,
                sequence: next_sequence,
            });
            next_sequence += 1;
        }
        subdivisions += 1;
        max_depth = max_depth.max(depth);
    }
    let upper = nodes.peek()?.upper;
    (lower <= upper && upper.is_finite()).then_some(MaximumBracket {
        lower,
        upper,
        subdivisions,
        max_depth,
    })
}
