//! Shared geometric proof of complete straight edge chains.
use crate::brep::{CurveSupport, OrientedEdgeUse};
use crate::math3::{cross, distance, dot, norm, normalize, sub};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;

const GEOM_TOL_MM: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIR_TOL: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

/// Each straight piece advances along one geometric segment. This is enough
/// to establish complete coverage of that segment, irrespective of splits.
pub fn line_chain(edges: &[&OrientedEdgeUse]) -> Option<f64> {
    let start = edges.first()?.start_mm;
    let end = edges.last()?.end_mm;
    if edges
        .windows(2)
        .any(|pair| distance(pair[0].end_mm, pair[1].start_mm) > GEOM_TOL_MM)
    {
        return None;
    }
    let direction = normalize(sub(end, start), DIR_TOL)?;
    let mut max_residual = 0.0_f64;
    for edge in edges {
        let CurveSupport::Line(line) = &edge.support else {
            return None;
        };
        let support_axis = normalize(line.direction, DIR_TOL)?;
        let move_axis = sub(edge.end_mm, edge.start_mm);
        if dot(move_axis, direction) <= GEOM_TOL_MM
            || dot(support_axis, direction).abs() < 1.0 - DIR_TOL
        {
            return None;
        }
        for point in [edge.start_mm, edge.end_mm] {
            let deviation = norm(cross(sub(point, start), direction))
                .max(norm(cross(sub(point, line.origin_mm), support_axis)));
            if deviation > GEOM_TOL_MM {
                return None;
            }
            max_residual = max_residual.max(deviation);
        }
    }
    Some(max_residual)
}
