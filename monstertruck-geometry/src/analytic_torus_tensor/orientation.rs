//! Certify the sign of the actual rational tensor's oriented normal on
//! every positive knot rectangle. Positional closeness alone cannot rule
//! out an arbitrarily narrow fold.
use super::ToroidalTensorChart;
use crate::{
    IntervalBezierSpan, IntervalControl, RationalControlBounds,
    analytic_bezier_partition::{interval_bezier_spans, same_algebraic_knot},
    prelude::*,
    proof_interval::Interval,
};

// Budgets bound proof work, not the geometry admitted by a successful proof.
const MAX_PATCHES: usize = 1 << 20;
const MAX_DEPTH: usize = 24;
type Net = Vec<Vec<IntervalControl>>;

fn clamped_continuous(knots: &[f64], degree: usize) -> bool {
    let Some((&first, &last)) = knots.first().zip(knots.last()) else {
        return false;
    };
    if knots.len() <= 2 * degree + 1
        || first >= last
        || !knots[..=degree]
            .iter()
            .all(|&x| same_algebraic_knot(x, first))
        || !knots[knots.len() - degree - 1..]
            .iter()
            .all(|&x| same_algebraic_knot(x, last))
    {
        return false;
    }
    let mut multiplicity = 1;
    for pair in knots[degree..knots.len() - degree].windows(2) {
        if same_algebraic_knot(pair[0], pair[1]) {
            multiplicity += 1;
            if multiplicity > degree {
                return false;
            }
        } else {
            multiplicity = 1;
        }
    }
    true
}

fn local_controls(net: &BsplineSurface<Vector4>, chart: ToroidalTensorChart) -> Option<Net> {
    let scale = net
        .control_points()
        .iter()
        .flatten()
        .map(|p| p.w)
        .reduce(f64::max)?;
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    let basis = [chart.radial, chart.axis.cross(chart.radial), chart.axis];
    let local = |p: &Vector4| -> Option<IntervalControl> {
        let weight = Interval::point(p.w).checked_div(Interval::point(scale))?;
        let xyz = [p.x, p.y, p.z];
        let offset: [Interval; 3] = std::array::from_fn(|i| {
            Interval::point(xyz[i]).minus(Interval::point(p.w).scale(chart.center[i]))
        });
        let coordinate = |i: usize| {
            (0..3)
                .map(|j| offset[j].scale(basis[i][j]))
                .fold(Interval::point(0.0), Interval::plus)
                .checked_div(Interval::point(scale))
        };
        Some([coordinate(0)?, coordinate(1)?, coordinate(2)?, weight])
    };
    net.control_points()
        .iter()
        .map(|row| row.iter().map(local).collect())
        .collect()
}

// Extract the tensor product of the same one-dimensional interval blossom
// used by the rational rim proof. Neither axis is sampled or knot-snapped.
fn spans_along_rows(net: &Net, knots: &[f64]) -> Option<Vec<Net>> {
    let mut spans = Vec::<Net>::new();
    for row in net {
        let extracted = interval_bezier_spans(row, knots)?;
        if spans.is_empty() {
            spans.resize_with(extracted.len(), Vec::new);
        }
        if spans.len() != extracted.len() {
            return None;
        }
        for (span, coefficients) in spans.iter_mut().zip(extracted) {
            span.push(coefficients.values);
        }
    }
    Some(spans)
}

fn transpose(net: &Net) -> Option<Net> {
    let count = net.first()?.len();
    if net.iter().any(|row| row.len() != count) {
        return None;
    }
    Some(
        (0..count)
            .map(|i| net.iter().map(|row| row[i]).collect())
            .collect(),
    )
}

fn split_rows(net: &Net) -> Option<[Net; 2]> {
    let mut first = Vec::with_capacity(net.len());
    let mut last = Vec::with_capacity(net.len());
    for row in net {
        let [a, b] = IntervalBezierSpan {
            values: row.clone(),
            parameter: [0.0, 1.0],
        }
        .split()?;
        first.push(a.values);
        last.push(b.values);
    }
    Some([first, last])
}

fn split(net: &Net, along_u: bool) -> Option<[Net; 2]> {
    if !along_u {
        return split_rows(net);
    }
    let [a, b] = split_rows(&transpose(net)?)?;
    Some([transpose(&a)?, transpose(&b)?])
}

fn oriented_normal(net: &Net, major: f64) -> Option<Interval> {
    let bounds = RationalControlBounds::from_controls(net.iter().flatten())?;
    let [x, y, z] = bounds.position;
    let rho = x.square().plus(y.square()).sqrt()?;
    if rho.lo <= 0.0 {
        return None;
    }
    let du = bounds.derivative(
        net.windows(2).flat_map(|rows| rows[0].iter().zip(&rows[1])),
        net.len() - 1,
    )?;
    let dv = bounds.derivative(
        net.iter()
            .flat_map(|row| row.windows(2).map(|p| (&p[0], &p[1]))),
        net.first()?.len() - 1,
    )?;
    // Outward torus normal field multiplied by the positive cylindrical
    // radius. Its dot with du × dv has the sign of the chart Jacobian.
    let radial = rho.minus(Interval::point(major));
    let normal = [x.times(radial), y.times(radial), z.times(rho)];
    let determinant = (0..3)
        .map(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            normal[i].times(du[j].times(dv[k]).minus(du[k].times(dv[j])))
        })
        .fold(Interval::point(0.0), Interval::plus);
    determinant.finite().then_some(determinant)
}

fn certify_patch(
    net: &Net,
    chart: ToroidalTensorChart,
    depth: usize,
    remaining: &mut usize,
) -> Option<()> {
    *remaining = remaining.checked_sub(1)?;
    if let Some(sign) = oriented_normal(net, chart.major_radius)
        && sign.directed()
    {
        return ((sign.lo > 0.0) == chart.outward).then_some(());
    }
    if depth == MAX_DEPTH {
        return None;
    }
    let [first, last] = split(net, depth.is_multiple_of(2))?;
    certify_patch(&first, chart, depth + 1, remaining)?;
    certify_patch(&last, chart, depth + 1, remaining)
}

pub(super) fn certify_orientation(
    net: &BsplineSurface<Vector4>,
    chart: ToroidalTensorChart,
) -> Option<()> {
    let u = net.knot_vector_u().as_ref();
    let v = net.knot_vector_v().as_ref();
    if !clamped_continuous(u, net.udegree()) || !clamped_continuous(v, net.vdegree()) {
        return None;
    }
    let local = local_controls(net, chart)?;
    let mut remaining = MAX_PATCHES;
    for u_span in spans_along_rows(&transpose(&local)?, u)? {
        for patch in spans_along_rows(&transpose(&u_span)?, v)? {
            certify_patch(&patch, chart, 0, &mut remaining)?;
        }
    }
    Some(())
}
