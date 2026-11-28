//! Remove artificial planar face partitions and redundant straight-edge vertices
//! from an already valid native shell. Used when a constructive selector spans
//! several exact Boolean-generated sections.
use crate::{Curve, Edge, EdgeId, Point3, Shell, Surface, VertexId, builder};
use monstertruck_core::tolerance::{DEFAULT_NUMERICAL_TOLERANCES, KERNEL_EPSILON};
use monstertruck_geometry::{
    analytic_carrier::CarrierDomainRelation, analytic_line_segment::SegmentLocus, prelude::*,
};
use monstertruck_topology::{
    compress::CompressedTrimmedShell,
    repair::{SurfaceDomainRelation, unify_same_domain_compressed_trimmed_faces},
    trimmed::TrimmedShell,
};
use rustc_hash::{FxHashMap as Map, FxHashSet as Set};

fn same_plane(a: &Surface, b: &Surface) -> Option<SurfaceDomainRelation> {
    if !matches!((a, b), (Surface::Plane(_), Surface::Plane(_))) {
        return None;
    }
    match crate::same_analytic_surface_domain(
        a,
        b,
        KERNEL_EPSILON,
        DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity,
    )? {
        CarrierDomainRelation::SameOrientation => Some(SurfaceDomainRelation::SameOrientation),
        CarrierDomainRelation::ReversedOrientation => {
            Some(SurfaceDomainRelation::ReversedOrientation)
        }
    }
}

fn certified_segment(edge: &Edge) -> bool {
    crate::certify_native_curve_segment_locus(
        &edge.curve(),
        edge.absolute_front().point(),
        edge.absolute_back().point(),
        KERNEL_EPSILON,
    ) == SegmentLocus::WithinTolerance
}

fn join_straight(a: &Edge, b: &Edge) -> Option<Edge> {
    let start = a.front().point();
    let middle = a.back().point();
    let end = b.back().point();
    let first = middle - start;
    let second = end - middle;
    let span = end - start;
    let length = span.norm();
    if !length.is_finite()
        || length <= KERNEL_EPSILON
        || first.dot(&second) <= 0.0
        || first.cross(&(span / length)).norm() > KERNEL_EPSILON
        || !certified_segment(a)
        || !certified_segment(b)
    {
        return None;
    }
    Some(builder::line(a.front(), b.back()))
}

fn collapse_line_vertices(shell: &mut Shell) {
    let mut incidence = Map::<VertexId, Vec<EdgeId>>::default();
    let mut seen = Set::default();
    for edge in shell.edge_iter() {
        if !seen.insert(edge.id()) {
            continue;
        }
        for vertex in [edge.front(), edge.back()] {
            incidence.entry(vertex.id()).or_default().push(edge.id());
        }
    }
    for (vertex, edges) in incidence {
        if edges.len() == 2 {
            shell.remove_vertex_by_join_edges(vertex, join_straight);
        }
    }
}

/// Normalize only certified coplanar face partitions and degree-two straight
/// seams. The complete result must remain a geometrically consistent closed shell.
/// # Errors
/// Returns a topology error if normalization cannot preserve the shell invariant.
pub fn normalize_planar_partitions(shell: &Shell) -> crate::Result<Shell> {
    let trimmed = TrimmedShell::<
        Point3,
        Curve,
        Surface,
        <Curve as ExactParameterBoundary2D<Surface>>::BoundaryCurve,
    >::from(shell.clone());
    let mut compressed = CompressedTrimmedShell::from(&trimmed);
    unify_same_domain_compressed_trimmed_faces(
        &mut compressed,
        same_plane,
        ExactParameterBoundary2D::exact_parameter_boundary_2d,
    )?;
    let trimmed = TrimmedShell::try_from(compressed)?;
    let mut normalized = trimmed.erase_trims();
    collapse_line_vertices(&mut normalized);
    normalized.check_solid_boundary()?;
    if normalized.is_geometric_consistent() {
        Ok(normalized)
    } else {
        Ok(shell.clone())
    }
}
