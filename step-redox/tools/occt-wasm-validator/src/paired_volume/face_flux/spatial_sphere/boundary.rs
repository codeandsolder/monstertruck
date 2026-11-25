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

/// The cumulative OCCT orientation identifies directed start/end vertices
/// even when a closed edge has two occurrences of the same vertex.
pub(super) fn directed_vertices(
    kernel: &mut OcctKernel,
    edge: ShapeHandle,
) -> Result<[ShapeHandle; 2]> {
    let mut vertices = [None; 2];
    for vertex in kernel.get_sub_shape_occurrences(edge, "vertex")? {
        let index = match kernel.shape_orientation(vertex)?.as_str() {
            "forward" => 0,
            "reversed" => 1,
            _ => bail!("non-oriented edge vertex"),
        };
        ensure!(
            vertices[index].replace(vertex).is_none(),
            "ambiguous directed edge vertex"
        );
    }
    let [start, end] = vertices;
    Ok([
        start.context("missing start vertex")?,
        end.context("missing end vertex")?,
    ])
}

pub(super) struct ClosedWire {
    pub(super) spans: Vec<Span>,
    pub(super) edges: usize,
    pub(super) closure_error: f64,
}

fn endpoints(span: &[Homogeneous]) -> Result<[Point; 2]> {
    Ok([
        cartesian(*span.first().context("missing span start")?),
        cartesian(*span.last().context("missing span end")?),
    ])
}

fn connect(spans: &mut Vec<Span>, start: Point, end: Point, error: &mut f64) -> Result<()> {
    let distance = norm(sub(end, start));
    ensure!(distance <= RESIDUAL_LIMIT, "open spatial boundary span");
    *error += distance;
    if distance > 0.0 {
        // Radially projecting this straight chord gives the short spherical
        // connector, including when the endpoints have tiny radial errors.
        spans.push(vec![homogeneous(start, 1.0), homogeneous(end, 1.0)]);
    }
    Ok(())
}

/// Find a complete directed endpoint matching rather than greedily
/// consuming the nearest end. A valid wire can revisit one vertex, in
/// which case several geometrically close assignments may be possible.
fn match_ends(options: &[Vec<usize>]) -> Result<Vec<usize>> {
    let count = options.len();
    let mut start_to_end = vec![None; count];
    let mut end_to_start = vec![None; count];
    for start in 0..count {
        // Breadth-first alternating paths avoid recursion on large wires.
        let mut queue = std::collections::VecDeque::from([start]);
        let mut seen_start = vec![false; count];
        seen_start[start] = true;
        let mut predecessor = vec![None; count];
        let free = 'search: loop {
            let current = queue
                .pop_front()
                .context("unmatched spatial boundary vertex")?;
            for &end in &options[current] {
                if predecessor[end].is_some() {
                    continue;
                }
                predecessor[end] = Some(current);
                if let Some(owner) = end_to_start[end] {
                    if !seen_start[owner] {
                        seen_start[owner] = true;
                        queue.push_back(owner);
                    }
                } else {
                    break 'search end;
                }
            }
        };
        let mut end = free;
        loop {
            let owner = predecessor[end].context("missing matching predecessor")?;
            let previous = start_to_end[owner].replace(end);
            end_to_start[end] = Some(owner);
            let Some(previous) = previous else {
                break;
            };
            end = previous;
        }
    }
    start_to_end
        .into_iter()
        .map(|end| end.context("unmatched directed start"))
        .collect()
}

/// Close the radially projected spatial wire explicitly, including between
/// consecutive Bézier spans. Native callers additionally require the same
/// topological vertex at every edge join.
pub(super) fn close_wire(
    edges: Vec<Vec<Span>>,
    mut same_vertex: impl FnMut(usize, usize) -> Result<bool>,
) -> Result<ClosedWire> {
    ensure!(!edges.is_empty(), "empty spherical wire");
    let ends = edges
        .iter()
        .map(|edge| -> Result<_> {
            Ok([
                endpoints(edge.first().context("missing first span")?)?[0],
                endpoints(edge.last().context("missing last span")?)?[1],
            ])
        })
        .collect::<Result<Vec<_>>>()?;
    let mut options = Vec::with_capacity(edges.len());
    for (start_index, [start, _]) in ends.iter().copied().enumerate() {
        let mut compatible = Vec::new();
        for (end_index, [_, end]) in ends.iter().copied().enumerate() {
            let distance = norm(sub(start, end));
            if distance <= RESIDUAL_LIMIT && same_vertex(end_index, start_index)? {
                compatible.push((end_index, distance));
            }
        }
        compatible.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        options.push(compatible.into_iter().map(|(index, _)| index).collect());
    }
    let matching = match_ends(&options)?;
    let edge_count = edges.len();
    let mut spans = Vec::new();
    let mut closure_error = 0.0;
    for edge in edges {
        let mut previous = None;
        for span in edge {
            let [start, end] = endpoints(&span)?;
            if let Some(previous) = previous {
                connect(&mut spans, previous, start, &mut closure_error)?;
            }
            spans.push(span);
            previous = Some(end);
        }
    }
    for (start_index, end_index) in matching.into_iter().enumerate() {
        connect(
            &mut spans,
            ends[end_index][1],
            ends[start_index][0],
            &mut closure_error,
        )?;
    }
    Ok(ClosedWire {
        spans,
        edges: edge_count,
        closure_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_finds_alternative_for_shared_vertex() -> Result<()> {
        assert_eq!(match_ends(&[vec![0, 1], vec![0]])?, [1, 0]);
        assert!(match_ends(&[vec![0], vec![0]]).is_err());
        Ok(())
    }

    #[test]
    fn closure_includes_internal_spans_and_respects_native_vertex_identity() -> Result<()> {
        let mut edge = circle_spans(
            [0.0; 3],
            1.0,
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            0.0,
            std::f64::consts::TAU,
        )?;
        edge[1][0][0] += 1e-13;
        let count = edge.len();
        let closed = close_wire(vec![edge.clone()], |_, _| Ok(true))?;
        assert!(closed.spans.len() > count && closed.closure_error > 0.0);
        assert!(close_wire(vec![edge], |_, _| Ok(false)).is_err());
        Ok(())
    }
}
