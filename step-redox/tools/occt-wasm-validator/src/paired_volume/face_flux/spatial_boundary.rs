//! Shared exact spatial edge representations and directed wire closure.
//! Each carrier integration independently checks its own incidence.

use super::{
    array,
    rational_bezier::{self, Homogeneous, spatial},
};
use crate::vector::{Point, add, cross, dot, norm, scale, sub, xyz};
use anyhow::{Context, Result, bail, ensure};
use occt_wasm::{CurveGeometry, OcctKernel, ShapeHandle};

const FRAME_TOLERANCE: f64 = 1e-11;
pub(super) type Span = Vec<Homogeneous>;

pub(super) fn cartesian(h: Homogeneous) -> Point { scale(spatial(h), h[3].recip()) }

pub(super) fn homogeneous(point: Point, weight: f64) -> Homogeneous {
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
    ellipse_spans(offset, [radius; 2], axis, x_direction, first, last)
}

fn ellipse_spans(
    offset: Point,
    radii: [f64; 2],
    axis: Point,
    x_direction: Point,
    first: f64,
    last: f64,
) -> Result<Vec<Span>> {
    ensure!(
        radii.iter().all(|r| r.is_finite() && *r > 0.0),
        "invalid elliptic edge radii"
    );
    ensure!(
        first.is_finite()
            && last.is_finite()
            && first < last
            && last - first <= std::f64::consts::TAU + FRAME_TOLERANCE,
        "invalid elliptic edge trim"
    );
    let axis = crate::vector::unit(axis, FRAME_TOLERANCE).context("degenerate conic axis")?;
    let x = crate::vector::unit(x_direction, FRAME_TOLERANCE)
        .context("degenerate conic x direction")?;
    ensure!(
        dot(axis, x).abs() <= FRAME_TOLERANCE,
        "nonorthogonal conic frame"
    );
    let y = cross(axis, x);
    let radial = |t: f64| {
        let (sin, cos) = t.sin_cos();
        add(scale(x, radii[0] * cos), scale(y, radii[1] * sin))
    };
    let point = |t| add(offset, radial(t));
    let mut spans = Vec::new();
    let mut start = first;
    while start < last {
        let end = (start + std::f64::consts::FRAC_PI_2).min(last);
        ensure!(
            end > start && spans.len() < 5,
            "invalid elliptic span partition"
        );
        let midpoint = f64::midpoint(start, end);
        let weight = ((end - start) * 0.5).cos();
        let middle = add(offset, scale(radial(midpoint), weight.recip()));
        spans.push(vec![
            homogeneous(point(start), 1.0),
            homogeneous(middle, weight),
            homogeneous(point(end), 1.0),
        ]);
        start = end;
    }
    Ok(spans)
}

/// Evaluate in a span-local homogeneous frame, retaining tiny curve
/// derivatives even when the spatial coordinates are much larger.
pub(super) struct SpanEvaluator {
    anchor: Point,
    controls: Vec<Homogeneous>,
    derivatives: Vec<Homogeneous>,
}

fn de_casteljau(controls: &[Homogeneous], parameter: f64) -> Homogeneous {
    // Normal STEP line, conic and low-order spline spans need no allocation
    // per quadrature node. Higher-degree curves retain the same algorithm.
    let mut local = [[0.0; 4]; 16];
    let mut heap;
    let work = if controls.len() <= local.len() {
        local[..controls.len()].copy_from_slice(controls);
        &mut local[..controls.len()]
    } else {
        heap = controls.to_vec();
        heap.as_mut_slice()
    };
    for count in (1..work.len()).rev() {
        for i in 0..count {
            let a = work[i];
            let b = work[i + 1];
            work[i] = std::array::from_fn(|j| (b[j] - a[j]).mul_add(parameter, a[j]));
        }
    }
    work[0]
}

impl SpanEvaluator {
    pub(super) fn new(span: &[Homogeneous]) -> Result<Self> {
        ensure!(
            span.len() >= 2
                && span
                    .iter()
                    .all(|h| h[3] > 0.0 && h.iter().all(|x| x.is_finite())),
            "invalid spatial span"
        );
        let anchor = cartesian(span[0]);
        let controls = span
            .iter()
            .map(|h| {
                [
                    (-anchor[0]).mul_add(h[3], h[0]),
                    (-anchor[1]).mul_add(h[3], h[1]),
                    (-anchor[2]).mul_add(h[3], h[2]),
                    h[3],
                ]
            })
            .collect::<Vec<_>>();
        let degree = (controls.len() - 1) as f64;
        // Differentiate the original homogeneous control net once. The
        // penultimate evaluation-level subtraction loses low-order bits on
        // short spans, particularly when the rational weights vary.
        let derivatives = controls
            .windows(2)
            .map(|pair| std::array::from_fn(|j| degree * (pair[1][j] - pair[0][j])))
            .collect::<Vec<_>>();
        ensure!(
            controls
                .iter()
                .chain(&derivatives)
                .flatten()
                .all(|x| x.is_finite()),
            "nonfinite local homogeneous span"
        );
        Ok(Self {
            anchor,
            controls,
            derivatives,
        })
    }

    pub(super) fn evaluate(&self, parameter: f64) -> (Point, Point) {
        let h = de_casteljau(&self.controls, parameter);
        let d = de_casteljau(&self.derivatives, parameter);
        let local = cartesian(h);
        let derivative = std::array::from_fn(|j| (-local[j]).mul_add(d[3], d[j]) / h[3]);
        (add(self.anchor, local), derivative)
    }
}

pub(super) fn edge_spans(
    kernel: &mut OcctKernel,
    edge: ShapeHandle,
    center: Point,
    radius: f64,
    tolerance: f64,
) -> Result<Vec<Span>> {
    let data = kernel.curve_data(edge)?;
    let (first, last) = (data.first, data.last);
    let mut spans = match data.geometry {
        CurveGeometry::Line { point, direction } => {
            let point = scale(sub(xyz(point), center), radius.recip());
            let direction = scale(xyz(direction), radius.recip());
            vec![vec![
                homogeneous(add(point, scale(direction, first)), 1.0),
                homogeneous(add(point, scale(direction, last)), 1.0),
            ]]
        }
        CurveGeometry::Ellipse {
            center: c,
            axis,
            major_direction,
            major_radius,
            minor_radius,
        } => ellipse_spans(
            scale(sub(xyz(c), center), radius.recip()),
            [major_radius / radius, minor_radius / radius],
            xyz(axis),
            xyz(major_direction),
            first,
            last,
        )?,
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
        norm(sub(cartesian(*begin), start)) <= tolerance
            && norm(sub(cartesian(*finish), end)) <= tolerance,
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

fn connect(
    spans: &mut Vec<Span>,
    start: Point,
    end: Point,
    error: &mut f64,
    tolerance: f64,
) -> Result<()> {
    let distance = norm(sub(end, start));
    ensure!(distance <= tolerance, "open spatial boundary span");
    *error += distance;
    if distance > 0.0 {
        // A spatial chord supplies a small explicit connector. Each carrier
        // projects it as part of its own boundary integral.
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

/// Close the spatial wire explicitly, including between
/// consecutive Bézier spans. Native callers additionally require the same
/// topological vertex at every edge join.
pub(super) fn close_wire(
    edges: Vec<Vec<Span>>,
    tolerance: f64,
    mut same_vertex: impl FnMut(usize, usize) -> Result<bool>,
) -> Result<ClosedWire> {
    ensure!(!edges.is_empty(), "empty spatial wire");
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
            if distance <= tolerance && same_vertex(end_index, start_index)? {
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
                connect(&mut spans, previous, start, &mut closure_error, tolerance)?;
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
            tolerance,
        )?;
    }
    Ok(ClosedWire {
        spans,
        edges: edge_count,
        closure_error,
    })
}

/// Connect each finite curve to its actual topological vertices, rather
/// than connecting two independently located curves directly. Opposite
/// coedges on adjacent faces then share the same spatial boundary chain.
fn anchor_edge(edge: &mut Vec<Span>, vertices: [Point; 2], tolerance: f64) -> Result<f64> {
    let start = endpoints(edge.first().context("missing first span")?)?[0];
    let end = endpoints(edge.last().context("missing last span")?)?[1];
    let mut before = Vec::new();
    let mut error = 0.0;
    connect(&mut before, vertices[0], start, &mut error, tolerance)?;
    before.append(edge);
    connect(&mut before, end, vertices[1], &mut error, tolerance)?;
    *edge = before;
    Ok(error)
}

#[derive(Clone, Copy)]
pub(super) enum Closure {
    Direct,
    NativeVertices,
}

/// Extract every directed native wire occurrence using topological identity
/// as well as endpoint proximity. No pcurve or meshing data enters this path.
pub(super) fn native_wires(
    kernel: &mut OcctKernel,
    face: ShapeHandle,
    center: Point,
    length: f64,
    tolerance: f64,
    closure: Closure,
) -> Result<Vec<ClosedWire>> {
    ensure!(
        length.is_finite() && length > 0.0 && tolerance.is_finite() && tolerance > 0.0,
        "invalid spatial normalization"
    );
    let mut wires = Vec::new();
    for wire in kernel.get_sub_shape_occurrences(face, "wire")? {
        let mut edges = Vec::new();
        let mut vertices = Vec::new();
        let mut anchor_error = 0.0;
        for edge in kernel.get_sub_shape_occurrences(wire, "edge")? {
            let ids = directed_vertices(kernel, edge)?;
            let mut spans = edge_spans(kernel, edge, center, length, tolerance)?;
            if matches!(closure, Closure::NativeVertices) {
                let [start, end] = ids.map(|id| -> Result<Point> {
                    Ok(scale(
                        sub(array(kernel.vertex_position(id)?, "native vertex")?, center),
                        length.recip(),
                    ))
                });
                anchor_error += anchor_edge(&mut spans, [start?, end?], tolerance)?;
            }
            vertices.push(ids);
            edges.push(spans);
        }
        let mut closed = close_wire(edges, tolerance, |end, start| {
            Ok(kernel.is_same(vertices[end][1], vertices[start][0])?)
        })?;
        closed.closure_error += anchor_error;
        wires.push(closed);
    }
    ensure!(!wires.is_empty(), "missing spatial boundary wires");
    Ok(wires)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_spans_preserve_derivatives_in_large_coordinate_frames() -> Result<()> {
        let base = 2.0_f64.powi(26);
        let step = 2.0_f64.powi(-20);
        let span = vec![
            homogeneous([base, 0., 0.], 1.),
            homogeneous([base + step, 0., 0.], 1.),
            homogeneous([base + 2. * step, 0., 0.], 1.),
            homogeneous([base + 3. * step, 0., 0.], 1.),
        ];
        let evaluator = SpanEvaluator::new(&span)?;
        for i in 0..=1000 {
            let t = f64::from(i) / 1000.;
            let (_, derivative) = evaluator.evaluate(t);
            assert_eq!(derivative[0], 3. * step);
        }
        // Positive rational weights alter parameter speed, not the line
        // locus. The local evaluation must preserve the small displacement
        // after dividing by the changing homogeneous denominator.
        let rational = vec![
            homogeneous([base, 0., 0.], 1.),
            homogeneous([base + step, 0., 0.], 0.5),
        ];
        let evaluator = SpanEvaluator::new(&rational)?;
        for i in 0..=1000 {
            let t = f64::from(i) / 1000.;
            let (_, derivative) = evaluator.evaluate(t);
            let expected = 0.5 * step / (-0.5_f64).mul_add(t, 1.).powi(2);
            assert!((derivative[0] - expected).abs() <= 8. * f64::EPSILON * expected);
        }
        Ok(())
    }

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
        let closed = close_wire(vec![edge.clone()], FRAME_TOLERANCE, |_, _| Ok(true))?;
        assert!(closed.spans.len() > count && closed.closure_error > 0.0);
        assert!(close_wire(vec![edge], FRAME_TOLERANCE, |_, _| Ok(false)).is_err());
        Ok(())
    }
}
