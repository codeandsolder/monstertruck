//! Exact spatial edge representations, normalized to a unit sphere. All
//! admitted edges must satisfy the whole-domain sphere incidence bound.

use super::super::{
    array,
    rational_bezier::{self, Homogeneous, spatial},
};
use crate::vector::{Point, add, cross, dot, norm, scale, sub, xyz};
use anyhow::{Context, Result, bail, ensure};
use occt_wasm::{CurveGeometry, OcctKernel, ShapeHandle};

pub(super) const RESIDUAL_LIMIT: f64 = 1e-11;
pub(super) type Span = Vec<Homogeneous>;

pub(super) fn cartesian(h: Homogeneous) -> Point { scale(spatial(h), h[3].recip()) }

fn homogeneous(point: Point, weight: f64) -> Homogeneous {
    [
        point[0] * weight,
        point[1] * weight,
        point[2] * weight,
        weight,
    ]
}

pub(super) fn circle_spans(
    offset: Point,
    radius: f64,
    axis: Point,
    x_direction: Point,
    first: f64,
    last: f64,
) -> Result<Vec<Span>> {
    ensure!(
        radius.is_finite() && radius > 0.0,
        "invalid circular edge radius"
    );
    ensure!(
        first.is_finite()
            && last.is_finite()
            && first < last
            && last - first <= std::f64::consts::TAU + RESIDUAL_LIMIT,
        "invalid circular edge trim"
    );
    let axis = crate::vector::unit(axis, RESIDUAL_LIMIT).context("degenerate circle axis")?;
    let x = crate::vector::unit(x_direction, RESIDUAL_LIMIT)
        .context("degenerate circle x direction")?;
    ensure!(
        dot(axis, x).abs() <= RESIDUAL_LIMIT,
        "nonorthogonal circle frame"
    );
    let y = cross(axis, x);
    let direction = |t: f64| {
        let (sin, cos) = t.sin_cos();
        add(scale(x, cos), scale(y, sin))
    };
    let point = |t| add(offset, scale(direction(t), radius));
    let mut spans = Vec::new();
    let mut start = first;
    while start < last {
        let end = (start + std::f64::consts::FRAC_PI_2).min(last);
        ensure!(
            end > start && spans.len() < 5,
            "invalid circular span partition"
        );
        let midpoint = f64::midpoint(start, end);
        let weight = ((end - start) * 0.5).cos();
        let middle = add(offset, scale(direction(midpoint), radius / weight));
        spans.push(vec![
            homogeneous(point(start), 1.0),
            homogeneous(middle, weight),
            homogeneous(point(end), 1.0),
        ]);
        start = end;
    }
    Ok(spans)
}

/// The homogeneous polynomial Q·Q-W² divided by W² is the normalized
/// squared-radius residual. Bernstein coefficient extrema and positive
/// weights bound the entire span, not just evaluated samples.
pub(super) fn radial_bound(span: &[Homogeneous]) -> Result<f64> {
    let residual = rational_bezier::normalized_squared_radius_bound(span)?;
    ensure!(
        residual.is_finite() && residual <= RESIDUAL_LIMIT,
        "edge not on sphere: normalized squared-radius residual {residual:e}"
    );
    Ok(residual)
}

pub(super) fn point_derivative(span: &[Homogeneous], parameter: f64) -> (Point, Point) {
    let mut poles = span.to_vec();
    let degree = (span.len() - 1) as f64;
    let interpolate = |a: Homogeneous, b: Homogeneous| {
        std::array::from_fn(|j| b[j].mul_add(parameter, a[j] * (1.0 - parameter)))
    };
    while poles.len() > 2 {
        for i in 0..poles.len() - 1 {
            poles[i] = interpolate(poles[i], poles[i + 1]);
        }
        poles.pop();
    }
    // The penultimate De Casteljau level gives the derivative as well as
    // the final point, with one work buffer and no separate derivative fit.
    let h = interpolate(poles[0], poles[1]);
    let d: Homogeneous = std::array::from_fn(|i| degree * (poles[1][i] - poles[0][i]));
    let point = cartesian(h);
    let derivative = scale(sub(spatial(d), scale(point, d[3])), h[3].recip());
    (point, derivative)
}

pub(super) fn edge_spans(
    kernel: &mut OcctKernel,
    edge: ShapeHandle,
    center: Point,
    radius: f64,
) -> Result<Vec<Span>> {
    let data = kernel.curve_data(edge)?;
    let (first, last) = (data.first, data.last);
    let mut spans = match data.geometry {
        CurveGeometry::Circle {
            center: c,
            radius: r,
            axis,
            x_direction,
        } => circle_spans(
            scale(sub(xyz(c), center), radius.recip()),
            r / radius,
            xyz(axis),
            xyz(x_direction),
            first,
            last,
        )?,
        CurveGeometry::BSpline(ref spline) => {
            rational_bezier::spans(spline, center, radius, first, last)?
        }
        _ => bail!("unsupported exact spatial curve carrier"),
    };
    // Independently compare the extracted curve's endpoints to the native
    // finite edge, including any shape location already applied by OCCT.
    let start = scale(
        sub(
            array(kernel.curve_point_at_param(edge, first)?, "edge start")?,
            center,
        ),
        radius.recip(),
    );
    let end = scale(
        sub(
            array(kernel.curve_point_at_param(edge, last)?, "edge end")?,
            center,
        ),
        radius.recip(),
    );
    let begin = spans
        .first()
        .and_then(|s| s.first())
        .context("missing span start")?;
    let finish = spans
        .last()
        .and_then(|s| s.last())
        .context("missing span end")?;
    ensure!(
        norm(sub(cartesian(*begin), start)) <= RESIDUAL_LIMIT
            && norm(sub(cartesian(*finish), end)) <= RESIDUAL_LIMIT,
        "extracted spatial curve endpoints do not match native edge"
    );
    match kernel.shape_orientation(edge)?.as_str() {
        "forward" => (),
        "reversed" => {
            spans.reverse();
            for span in &mut spans {
                span.reverse();
            }
        }
        _ => bail!("non-oriented spatial edge"),
    }
    Ok(spans)
}

/// A directed wire's end multiset must match its start multiset. The
/// integral is additive and does not depend on OCCT's occurrence order.
pub(super) fn check_closed(edges: &[Vec<Span>]) -> Result<f64> {
    let endpoints = edges
        .iter()
        .map(|edge| -> Result<_> {
            let start = edge
                .first()
                .and_then(|s| s.first())
                .context("missing start")?;
            let end = edge.last().and_then(|s| s.last()).context("missing end")?;
            Ok((cartesian(*start), cartesian(*end)))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut ends = endpoints.iter().map(|&(_, end)| end).collect::<Vec<_>>();
    let mut closure_error = 0.0;
    for (start, _) in endpoints {
        let index = ends
            .iter()
            .position(|&end| norm(sub(start, end)) <= RESIDUAL_LIMIT)
            .context("open spatial boundary wire")?;
        closure_error += norm(sub(start, ends.swap_remove(index)));
    }
    ensure!(ends.is_empty(), "unmatched wire endpoints");
    Ok(closure_error)
}
