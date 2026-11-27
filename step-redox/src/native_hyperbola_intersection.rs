//! Native STEP INTERSECTION_CURVE proposal for a separately healed edge.
//!
//! The 3D master curve is an exact analytic HYPERBOLA; associated geometry
//! consists of its corrected plane and independently authored cone support.
//! The source spline and all original topology are retained elsewhere.
//! In particular the nonlinear cone-UV function is never exported as a
//! falsely exact rational PCURVE.

use crate::source_edge_healing::HyperbolaHealingCandidate;
use monstertruck_geometry::prelude::{
    EuclideanSpace, Line, Matrix4, Plane, Point3, Processor, RevolutionSurface, Vector3,
};
use monstertruck_io::step::load::step_geometry::{
    Conic3D, Curve3D, ElementarySurface, Surface, SurfaceCurve3D, SurfaceCurveAssociatedGeometry,
    SurfaceCurveKind, SurfaceCurveRepresentation,
};

/// Produce an exact, self-contained analytic STEP intersection *geometry*.
///
/// This has no `BRep` topological edges or faces. The corrected plane is
/// intentionally not claimed to be the identical manufacturer source plane:
/// its whole-trim mismatch to that authored plane was bounded when the
/// healing candidate was constructed. The source cone's reference direction
/// is preserved by the verified analytical UV chart.
///
/// Both supporting surfaces are emitted as SURFACE associations, not as
/// fabricated STEP PCURVEs. Native STEP consumers may choose a different
/// surface UV branch, but the 3D geometry remains the preferred master.
#[must_use = "The constructed exact, self-contained analytic STEP intersection *geometry* is returned rather than stored."]
pub fn native_healed_intersection_curve(candidate: &HyperbolaHealingCandidate) -> Curve3D {
    let analytic = candidate.analytic;
    let axis = analytic.cone_axis();
    let normal = analytic.plane_normal();
    let transverse = axis.cross(normal);
    let origin = analytic.reference_origin();
    let point_on_plane = origin + normal * analytic.signed_plane_offset();
    let plane = Plane::new(
        point_on_plane,
        point_on_plane + transverse,
        point_on_plane + axis,
    );
    let (source_x, source_y) = candidate.cone_uv.reference_axes();
    let frame = Matrix4::from_cols(
        source_x.extend(0.0),
        source_y.extend(0.0),
        axis.extend(0.0),
        origin.to_vec().extend(1.0),
    );
    let r0 = analytic.reference_radius();
    let slope = analytic.slope();
    let generator = Line(Point3::new(r0, 0.0, 0.0), Point3::new(r0 + slope, 0.0, 1.0));
    let cone = Processor::with_transform(
        RevolutionSurface::by_revolution(generator, Point3::new(0.0, 0.0, 0.0), Vector3::unit_z()),
        frame,
    );

    let leader = Curve3D::Conic(Conic3D::Hyperbola(candidate.modeling_curve()));
    Curve3D::SurfaceCurve(SurfaceCurve3D::new(
        SurfaceCurveKind::Intersection,
        Box::new(leader),
        vec![
            SurfaceCurveAssociatedGeometry::Surface(Box::new(Surface::Elementary(
                ElementarySurface::Plane(plane),
            ))),
            SurfaceCurveAssociatedGeometry::Surface(Box::new(Surface::Elementary(
                ElementarySurface::ConicalSurface(cone),
            ))),
        ],
        SurfaceCurveRepresentation::Curve3D,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        source_edge_healing::{HyperbolaHealingRequest, propose_plane_cone_hyperbola_healing},
        source_edge_locus::StepEdge,
    };
    use anyhow::{Result, bail};
    use monstertruck_geometry::prelude::{InnerSpace, ParametricCurve, ParametricSurface};
    use monstertruck_io::step::{
        load::{SurfaceCurveHolder, Table, step_geometry::SurfaceCurveAssociatedGeometry},
        save::StepDisplay,
    };
    use step_p21::tables::EntityTable;

    fn original_manufacturer() -> Result<Vec<HyperbolaHealingCandidate>> {
        let text = std::str::from_utf8(include_bytes!(
            "../validation/fixtures/production_radial_slot.step"
        ))?;
        let exchange = step_p21::parser::parse(text)?;
        let [section] = exchange.data.as_slice() else {
            bail!("manufacturer input must have exactly one DATA section");
        };
        [(93481, 125603, 58899), (120654, 75450, 36791)]
            .into_iter()
            .map(|(edge, plane, cone)| {
                propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
                    solid_id: 15385,
                    source: StepEdge {
                        id: edge,
                        entities: &section.entities,
                    },
                    source_plane_face_id: plane,
                    source_cone_face_id: cone,
                    source_coordinate_scale_mm: 1.0,
                    endpoint_tolerance_mm: 1.0e-8,
                    target_gap_mm: 1.0e-9,
                })
            })
            .collect()
    }

    #[test]
    fn both_manufacturer_edges_emit_independent_exact_plane_and_cone_supports() -> Result<()> {
        for candidate in original_manufacturer()? {
            let geometry = native_healed_intersection_curve(&candidate);
            let Curve3D::SurfaceCurve(intersection) = &geometry else {
                bail!("native healed curve must retain associated geometry");
            };
            assert_eq!(intersection.kind(), SurfaceCurveKind::Intersection);
            assert_eq!(
                intersection.master_representation(),
                SurfaceCurveRepresentation::Curve3D
            );
            let [
                SurfaceCurveAssociatedGeometry::Surface(plane),
                SurfaceCurveAssociatedGeometry::Surface(cone),
            ] = intersection.associated_geometry()
            else {
                bail!("native healed curve requires exactly two analytic surfaces");
            };
            let Surface::Elementary(ElementarySurface::Plane(plane)) = plane.as_ref() else {
                bail!("first support must be the corrected analytic plane");
            };
            let Surface::Elementary(ElementarySurface::ConicalSurface(cone)) = cone.as_ref() else {
                bail!("second support must be the original-cone analytic geometry");
            };
            // Production computations are derived from exact formulae.
            // These evaluations independently guard against frame/axis or
            // trim changes in the native STEP emission code.
            for i in 0..=256 {
                let t = f64::from(i).mul_add(
                    (candidate.trim.end_parameter - candidate.trim.start_parameter) / 256.0,
                    candidate.trim.start_parameter,
                );
                let point = candidate
                    .analytic
                    .evaluate(t)
                    .ok_or_else(|| anyhow::anyhow!("invalid native hyperbola parameter"))?;
                let uv = candidate
                    .cone_uv
                    .evaluate(t)
                    .ok_or_else(|| anyhow::anyhow!("invalid cone UV parameter"))?;
                let mapped = cone.evaluate(uv.x, uv.y);
                assert!((mapped - point).magnitude() < 1.0e-10);
                let residual = (point - plane.origin()).dot(plane.normal()).abs();
                assert!(residual < 1.0e-10);
            }
            let fragment = StepDisplay::new(&geometry, 1).to_string();
            assert!(fragment.starts_with("#1 = INTERSECTION_CURVE("));
            assert!(fragment.contains("HYPERBOLA("));
            assert!(fragment.contains("PLANE("));
            assert!(fragment.contains("CONICAL_SURFACE("));
            assert!(fragment.contains(".CURVE_3D."));
            assert!(!fragment.contains("PCURVE("));
            let full = format!(
                "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('analytic healing'),'2;1');\nFILE_NAME('healed','2026-11-18',(''),(''),'','', '');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{fragment}\nENDSEC;\nEND-ISO-10303-21;\n"
            );
            let parsed = step_p21::parser::parse(&full)?;
            assert_eq!(parsed.data.len(), 1);
            let table = Table::from_step(&full)?;
            assert_eq!(table.surface_curve.len(), 1);
            assert_eq!(table.hyperbola.len(), 1);
            assert_eq!(table.plane.len(), 1);
            assert_eq!(table.conical_surface.len(), 1);
            assert!(
                table.entity_report.is_empty(),
                "native loader must parse every geometric record"
            );
            // Resolve the actual typed STEP entity graph rather than just
            // counting syntactically recognized records. Two SURFACE
            // associations must reconstruct an IntersectionCurve, not a
            // lossy SurfaceCurve or an invented PCURVE.
            let owned = <Table as EntityTable<SurfaceCurveHolder>>::get_owned(&table, 1)?;
            let loaded = Curve3D::try_from(&owned)
                .map_err(|error| anyhow::anyhow!("native intersection import: {error}"))?;
            let Curve3D::IntersectionCurve(imported) = loaded else {
                bail!("two native STEP surface associations lost their intersection type");
            };
            assert!(matches!(
                imported.leader().as_ref(),
                Curve3D::Conic(Conic3D::Hyperbola(_))
            ));
            assert!(matches!(
                imported.surface0().as_ref(),
                Surface::Elementary(ElementarySurface::Plane(_))
            ));
            assert!(matches!(
                imported.surface1().as_ref(),
                Surface::Elementary(ElementarySurface::ConicalSurface(_))
            ));
            // STEP's HYPERBOLA entity has no finite trim of its own, but
            // the imported native master must retain the same parameter
            // meaning and lie on both associated supports. The original
            // EDGE_CURVE vertices will define the trim in a later BRep.
            for i in 0..=8 {
                let parameter = f64::from(i).mul_add(
                    (candidate.trim.end_parameter - candidate.trim.start_parameter) / 8.0,
                    candidate.trim.start_parameter,
                );
                let original = candidate.analytic.evaluate(parameter).ok_or_else(|| {
                    anyhow::anyhow!("invalid source analytic hyperbola parameter")
                })?;
                let reimported = imported.evaluate(parameter);
                assert!(
                    (reimported - original).magnitude() < 1.0e-9,
                    "typed native re-import must retain the exact hyperbola and both carriers"
                );
            }
        }
        Ok(())
    }
}
