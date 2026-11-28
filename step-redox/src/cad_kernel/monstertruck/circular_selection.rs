//! Bounded circular selector: whole-curve support and directed arc coverage.
use super::{EDGE_SELECTOR_TOLERANCE_MM, segment_point};
use anyhow::{Context, Result, ensure};
use monstertruck_geometry::{
    analytic_curve_incidence::{CircularArcFrame, certify_circular_arc_coverage},
    prelude::TryIntoHomogeneousBsplineCurve,
};
use monstertruck_modeling::{Edge, InnerProduct, NurbsCurve, Point3, Shell, Vector3};
use std::collections::{HashMap, HashSet};

const TURN: f64 = std::f64::consts::TAU;

/// Prove a single connected, non-branching open chain using topology IDs,
/// not a coordinate weld. Check that its terminal points are the selector's
/// requested endpoints.
fn verify_chain(edges: &[Edge], frame: CircularArcFrame, sweep: f64) -> Result<()> {
    let start = frame.point(0.0);
    let end = frame.point(sweep);
    let mut adjacency = HashMap::<_, Vec<_>>::new();
    let mut points = HashMap::new();
    for edge in edges {
        let a = edge.absolute_front();
        let b = edge.absolute_back();
        ensure!(
            a.id() != b.id(),
            "bounded circular chain contains a closed edge"
        );
        points.insert(a.id(), a.point());
        points.insert(b.id(), b.point());
        adjacency.entry(a.id()).or_default().push(b.id());
        adjacency.entry(b.id()).or_default().push(a.id());
    }
    let endpoints = adjacency
        .iter()
        .filter(|(_, neighbors)| neighbors.len() == 1)
        .map(|(id, _)| *id)
        .collect::<Vec<_>>();
    let [a, b] = endpoints.as_slice() else {
        anyhow::bail!("bounded circular selector does not resolve to one open chain");
    };
    ensure!(
        adjacency
            .values()
            .all(|neighbors| (1..=2).contains(&neighbors.len())),
        "bounded circular selector has a branching chain"
    );
    let p = points[a];
    let q = points[b];
    let close = |p: Point3, q: Point3| (p - q).norm() <= EDGE_SELECTOR_TOLERANCE_MM;
    ensure!(
        (close(p, start) && close(q, end)) || (close(q, start) && close(p, end)),
        "bounded circular chain has incorrect endpoints"
    );
    let mut seen = HashSet::new();
    let mut stack = vec![*a];
    while let Some(v) = stack.pop() {
        if seen.insert(v) {
            stack.extend(adjacency[&v].iter().copied());
        }
    }
    ensure!(
        seen.len() == adjacency.len(),
        "bounded circular selector contains disconnected edges"
    );
    Ok(())
}

pub fn select_arc(
    shell: &Shell,
    center: [f64; 3],
    normal: [f64; 3],
    radius: f64,
    start: [f64; 3],
    sweep: f64,
) -> Result<Vec<Edge>> {
    let frame = CircularArcFrame::new(
        segment_point(center),
        Vector3::new(normal[0], normal[1], normal[2]),
        radius,
        segment_point(start),
        EDGE_SELECTOR_TOLERANCE_MM,
    )
    .context("invalid circular selector frame")?;
    let tolerance = frame.angle_tolerance();
    ensure!(
        sweep.is_finite() && sweep > tolerance && sweep < TURN - tolerance,
        "invalid bounded circular selector sweep"
    );
    let mut unique = HashSet::new();
    let mut matched = Vec::new();
    for edge in shell.edge_iter() {
        if !unique.insert(edge.id()) {
            continue;
        }
        let Some(coverage) = edge
            .oriented_curve()
            .try_into_homogeneous_bspline_curve()
            .and_then(|curve| certify_circular_arc_coverage(&NurbsCurve::new(curve), frame))
        else {
            continue;
        };
        if coverage.start_rad + coverage.sweep_rad <= sweep + tolerance {
            matched.push((coverage, edge));
        }
    }
    matched.sort_by(|(a, _), (b, _)| a.start_rad.total_cmp(&b.start_rad));
    let mut cursor = 0.0;
    for (coverage, _) in &matched {
        ensure!(
            (coverage.start_rad - cursor).abs() <= tolerance,
            "bounded circular selector coverage has a gap or overlap"
        );
        cursor = coverage.start_rad + coverage.sweep_rad;
    }
    ensure!(
        (cursor - sweep).abs() <= tolerance,
        "bounded circular selector does not cover the requested arc"
    );
    let edges = matched
        .into_iter()
        .map(|(_, edge)| edge)
        .collect::<Vec<_>>();
    verify_chain(&edges, frame, sweep).context("validate bounded circular selector topology")?;
    Ok(edges)
}
