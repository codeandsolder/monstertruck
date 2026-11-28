//! Three-dimensional incidence certificates for every boundary of a rewritten
//! STEP plane face. The source representation uncertainty is not enlarged.
use crate::{
    brep::{self, CurveSupport, OrientedEdgeUse, PlaneSupport},
    step_graph::{entity_ref_value, simple_record},
};
use anyhow::{Context, Result, bail, ensure};
use monstertruck_geometry::{
    analytic_plane_deviation::bound_nurbs_plane_deviation,
    prelude::{
        BoundedCurve, BsplineCurve, InnerProduct, KnotVector, Line, Matrix4, NurbsCurve, Point3,
        TryIntoHomogeneousBsplineCurve, Vector3, Vector4,
    },
};
use monstertruck_io::step::load::{
    EdgeCurveHolder, Table,
    step_geometry::{Conic3D, Curve3D},
};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use step_p21::{
    ast::{DataSection, EntityInstance, Parameter},
    tables::EntityTable,
};

/// Complete 3D boundary-incidence evidence, not a whole-solid validity or
/// parameter-space trimming certificate.
#[derive(Clone, Debug, Serialize)]
pub struct PlaneBoundaryEvidence {
    pub loop_count: usize,
    pub retained_edge_count: usize,
    pub retained_max_residual_upper_mm: f64,
    pub vertex_max_residual_upper_mm: f64,
    pub healed_edge_residual_upper_mm: f64,
    pub source_tolerance_mm: f64,
}

#[derive(Clone, Copy)]
struct PlaneBound {
    origin: Point3,
    normal: Vector3,
    gap: f64,
    tolerance: f64,
}

const fn point(p: [f64; 3]) -> Point3 { Point3::new(p[0], p[1], p[2]) }

impl PlaneBound {
    fn spline(self, curve: BsplineCurve<Vector4>) -> Result<f64> {
        let bracket = bound_nurbs_plane_deviation(
            &NurbsCurve::new(curve),
            self.origin,
            self.normal,
            self.gap,
        )
        .context("cannot bound entire retained curve against replacement plane")?;
        ensure!(
            bracket.upper <= self.tolerance,
            "boundary plane incidence not certified: upper deviation {} exceeds source budget {} (source units)",
            bracket.upper,
            self.tolerance
        );
        Ok(bracket.upper)
    }

    fn segment(self, a: Point3, b: Point3) -> Result<f64> {
        self.spline(
            Line(a, b)
                .try_into_homogeneous_bspline_curve()
                .context("cannot construct exact line segment")?,
        )
    }

    // Bound the *authored* line interval, not the importer's synthesized chord
    // through two approximate STEP vertices. Include both projected endpoints.
    fn line(self, edge: &OrientedEdgeUse, line: brep::LineSupport) -> Result<f64> {
        let origin = point(line.origin_mm);
        let direction = Vector3::new(line.direction[0], line.direction[1], line.direction[2]);
        let project = |p: Point3| -> Result<(Point3, f64)> {
            let projected = origin + direction * (p - origin).dot(&direction);
            // The line placement can be far from the finite edge. Account
            // for cancellation in the projection, not only the final chord.
            let coordinate_scale = [origin.x, origin.y, origin.z, p.x, p.y, p.z]
                .into_iter()
                .map(f64::abs)
                .fold(1.0_f64, f64::max);
            let margin = 256.0 * f64::EPSILON * coordinate_scale;
            ensure!(
                (projected - p).norm() + margin <= self.tolerance,
                "source vertex is not certified on its authored line within source uncertainty"
            );
            Ok((projected, margin))
        };
        let (start, start_margin) = project(point(edge.start_mm))?;
        let (end, end_margin) = project(point(edge.end_mm))?;
        let upper = self.segment(start, end)? + start_margin.max(end_margin);
        ensure!(
            upper <= self.tolerance,
            "projected source line exceeds replacement plane budget"
        );
        Ok(upper)
    }
}

// A rational quadratic exactly represents a finite conic segment. For a
// hyperbola, the homogeneous middle control is the mid-parameter point with
// weight cosh(half_span); the endpoints have unit weight. A parabola is the
// ordinary polynomial quadratic. These controls cover the entire trim.
fn conic_net(
    transform: Matrix4,
    range: (f64, f64),
    hyperbola: bool,
) -> Result<BsplineCurve<Vector4>> {
    let (a, b) = range;
    ensure!(
        a.is_finite() && b.is_finite() && b > a,
        "invalid conic trim"
    );
    let controls = if hyperbola {
        let mid = a.midpoint(b);
        [
            Vector4::new(a.cosh(), a.sinh(), 0.0, 1.0),
            Vector4::new(mid.cosh(), mid.sinh(), 0.0, ((b - a) * 0.5).cosh()),
            Vector4::new(b.cosh(), b.sinh(), 0.0, 1.0),
        ]
    } else {
        [
            Vector4::new(a * a, 2.0 * a, 0.0, 1.0),
            Vector4::new(a * b, a + b, 0.0, 1.0),
            Vector4::new(b * b, 2.0 * b, 0.0, 1.0),
        ]
    };
    Ok(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(2),
        controls.into_iter().map(|p| transform * p).collect(),
    ))
}

fn curve_upper(curve: &Curve3D, plane: PlaneBound) -> Result<f64> {
    // A declared 3D master is validated from the original entity graph
    // before this function is invoked. Its leader, not a projected kernel
    // approximation, is the authored curve.
    match curve {
        Curve3D::SurfaceCurve(curve) => curve_upper(curve.leader(), plane),
        Curve3D::IntersectionCurve(curve) => curve_upper(curve.leader(), plane),
        Curve3D::Line(_) => bail!("retained authored line lacks a valid source support"),
        Curve3D::Polyline(polyline) => {
            ensure!(polyline.len() >= 2, "empty or degenerate source polyline");
            polyline.windows(2).try_fold(0.0_f64, |maximum, points| {
                Ok(maximum.max(plane.segment(points[0], points[1])?))
            })
        }
        Curve3D::Conic(Conic3D::Hyperbola(conic)) => {
            plane.spline(conic_net(*conic.transform(), conic.range_tuple(), true)?)
        }
        Curve3D::Conic(Conic3D::Parabola(conic)) => {
            plane.spline(conic_net(*conic.transform(), conic.range_tuple(), false)?)
        }
        _ => plane.spline(
            curve
                .try_into_homogeneous_bspline_curve()
                .context("no exact curve-locus representation for retained source boundary")?,
        ),
    }
}

// Surface/ intersection/seam curve wrappers can declare a PCURVE, not the 3D
// leader, as their master. A leader-only certificate cannot authorize those.
fn require_3d_master(
    mut id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Result<u64> {
    let mut seen = HashSet::new();
    loop {
        ensure!(seen.insert(id), "cyclic source curve master references");
        let entity = index
            .get(&id)
            .and_then(|&i| entities.get(i))
            .context("missing source boundary curve")?;
        let Some(record) = simple_record(entity) else {
            return Ok(id);
        };
        if !matches!(
            record.name.as_str(),
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE"
        ) {
            ensure!(
                record.name != "PCURVE",
                "source boundary requires a parameter-space master proof"
            );
            return Ok(id);
        }
        let Parameter::List(params) = &record.parameter else {
            bail!("malformed source surface curve");
        };
        ensure!(
            matches!(params.get(3), Some(Parameter::Enumeration(master)) if master == "CURVE_3D"),
            "source boundary does not declare a 3D master curve"
        );
        id = params
            .get(1)
            .and_then(entity_ref_value)
            .context("source surface curve lacks a 3D leader")?;
    }
}

fn retained_edge_upper(
    source: &DataSection,
    index: &HashMap<u64, usize>,
    table: &Table,
    edge: &OrientedEdgeUse,
    plane: PlaneBound,
) -> Result<f64> {
    let leader = require_3d_master(edge.curve_id, &source.entities, index)?;
    if let CurveSupport::Line(line) = brep::curve_support(leader, &source.entities, index) {
        return plane.line(edge, line);
    }
    let owned = <Table as EntityTable<EdgeCurveHolder>>::get_owned(table, edge.edge_id)?;
    let curve = owned
        .parse_curve3d()
        .map_err(|error| anyhow::anyhow!("{error}"))?;
    curve_upper(&curve, plane)
}

pub(super) fn certify(
    source: &DataSection,
    index: &HashMap<u64, usize>,
    request: super::FullSourceHyperbolaRepair,
    candidate: &crate::source_edge_healing::HyperbolaHealingCandidate,
    replacement: PlaneSupport,
    scale: f64,
) -> Result<PlaneBoundaryEvidence> {
    let source_tolerance_mm = candidate.evidence.source_representation_tolerance_mm;
    let plane = PlaneBound {
        origin: point(replacement.origin_mm),
        normal: Vector3::new(
            replacement.normal[0],
            replacement.normal[1],
            replacement.normal[2],
        ),
        gap: request.target_gap_mm / scale,
        tolerance: source_tolerance_mm / scale,
    };
    ensure!(
        [plane.gap, plane.tolerance]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0),
        "source units cannot represent boundary incidence budget"
    );
    let loops = brep::face_loops(request.plane_face_id, &source.entities, index)
        .context("cannot enumerate every replacement plane boundary")?;
    ensure!(!loops.is_empty(), "replacement plane has no boundary loops");
    let table = Table::from_data_section(source);
    let mut checked = HashSet::new();
    let mut retained_max = 0.0_f64;
    let mut vertex_max = 0.0_f64;
    for boundary in &loops {
        ensure!(
            !boundary.edges.is_empty(),
            "replacement plane has an empty boundary loop"
        );
        for edge in &boundary.edges {
            // The affine segment joining the two vertices is in the plane's
            // slab iff both points are; the same global bound supplies the
            // numerical allowance as the other curve certificates.
            vertex_max = vertex_max.max(
                plane
                    .segment(point(edge.start_mm), point(edge.end_mm))
                    .with_context(|| {
                        format!(
                            "replacement face #{} edge #{} vertices",
                            request.plane_face_id, edge.edge_id
                        )
                    })?,
            );
            if edge.edge_id != request.edge_id && checked.insert(edge.edge_id) {
                retained_max = retained_max.max(
                    retained_edge_upper(source, index, &table, edge, plane).with_context(|| {
                        format!(
                            "replacement face #{} retained edge #{}",
                            request.plane_face_id, edge.edge_id
                        )
                    })?,
                );
            }
        }
    }
    let healed_edge_residual_upper_mm = candidate
        .analytic
        .source_plane_residual_upper(plane.origin * scale, plane.normal, candidate.trim)
        .context("cannot bound healed hyperbola on serialized replacement plane")?;
    ensure!(
        healed_edge_residual_upper_mm <= source_tolerance_mm,
        "healed hyperbola exceeds replacement plane source incidence budget"
    );
    Ok(PlaneBoundaryEvidence {
        loop_count: loops.len(),
        retained_edge_count: checked.len(),
        retained_max_residual_upper_mm: retained_max * scale,
        vertex_max_residual_upper_mm: vertex_max * scale,
        healed_edge_residual_upper_mm,
        source_tolerance_mm,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entire_boundary_rejects_remote_vertices_and_interior_spline_excursion() {
        let plane = PlaneBound {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::z(),
            tolerance: 1.0e-5,
            gap: 1.0e-9,
        };
        assert!(
            plane
                .segment(Point3::new(0.0, 0.0, 0.0), Point3::new(100.0, 0.0, 2.0e-5))
                .is_err()
        );
        let curve = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 1.0),
                Vector4::new(0.5, 0.0, 1.0e-3, 1.0),
                Vector4::new(1.0, 0.0, 0.0, 1.0),
            ],
        );
        assert!(plane.spline(curve).is_err());
    }

    #[test]
    fn exact_conic_nets_include_entire_locus() -> Result<()> {
        use monstertruck_geometry::prelude::ParametricCurve;
        let (a, b) = (-0.7_f64, 1.3_f64);
        let transform = Matrix4::new_scaling(2.0);
        for hyperbola in [false, true] {
            let net = NurbsCurve::new(conic_net(transform, (a, b), hyperbola)?);
            for i in 0..=20 {
                let p = net.evaluate(f64::from(i) / 20.0);
                let residual = if hyperbola {
                    (p.x * p.x - p.y * p.y - 4.0).abs()
                } else {
                    (p.y * p.y - 8.0 * p.x).abs()
                };
                assert!(residual < 1.0e-12, "conic residual {residual}");
            }
        }
        Ok(())
    }
}
