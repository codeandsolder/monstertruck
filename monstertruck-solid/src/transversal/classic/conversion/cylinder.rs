//! Replace a regular marching intersection with a bounded analytic leader.
//!
//! The marching polyline selects the finite angular interval and transverse
//! root, not the exported geometry. Both exact cylinder carriers must be
//! recognized independently; the whole polyline must cover the chosen branch
//! within the mesh discretization budget. The cubic leader itself is generated
//! to the separate structural precision budget.

use super::{BooleanTolerance, PolylineCurve, polyline_leader_bspline, prove_iso_candidate};
use crate::transversal::classic::ShapeOpsSurface;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use monstertruck_geometry::{
    analytic_cylinder_intersection::{BoundedSourceIntersection, CylinderCarrier},
    prelude::*,
};

fn carrier<S: ShapeOpsSurface>(surface: &S) -> Option<CylinderCarrier> {
    CylinderCarrier::from_homogeneous_surface(surface.try_into_homogeneous_bspline_surface()?)
}

/// Construct an independently bounded leader in the original traversal
/// direction. Finite chart checks are additional guards: the analytic
/// incidence and branch coverage do not depend on sampled validation.
pub(super) fn crossing_cylinder_leader<S: ShapeOpsSurface>(
    intersection: &IntersectionCurve<PolylineCurve<Point3>, S, S>,
    tolerances: BooleanTolerance,
) -> Option<NurbsCurve<Vector4>> {
    let first = carrier(intersection.surface0())?;
    let second = carrier(intersection.surface1())?;
    let guide = intersection.leader();
    let source = NurbsCurve::from(polyline_leader_bspline(guide)?);
    let precision = tolerances
        .length
        .min(DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_structural_length);
    let mut selected = None::<BoundedSourceIntersection>;
    for (primary, secondary) in [(first, second), (second, first)] {
        let Some(bound) = BoundedSourceIntersection::from_spline(
            primary,
            secondary,
            &source,
            precision * 0.5,
            tolerances.mesh_seed / 16.0,
        ) else {
            continue;
        };
        if bound.separation_upper <= tolerances.mesh_seed
            && selected.is_none_or(|old| bound.separation_upper < old.separation_upper)
        {
            selected = Some(bound);
        }
    }
    let arc = selected?.proposal.arc;
    let leader = arc.cubic_leader(
        precision * 0.5,
        DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_leader_max_segments,
    )?;
    let mut curve = leader.curve;
    let endpoint_error = curve
        .front()
        .distance(guide.front())
        .max(curve.back().distance(guide.back()));
    if !endpoint_error.is_finite() || leader.maximum_error + endpoint_error > precision {
        return None;
    }
    // A clamped polynomial spline has nonnegative partition-of-unity basis
    // functions. Moving only the two endpoint controls changes the whole
    // locus by at most their maximum displacement, already budgeted above.
    let last = curve.control_points().len().checked_sub(1)?;
    *curve.control_point_mut(0) = guide.front();
    *curve.control_point_mut(last) = guide.back();
    let result = NurbsCurve::from(curve);
    prove_iso_candidate(
        &result,
        intersection.surface0(),
        intersection.surface1(),
        guide,
        tolerances,
    )
    .then_some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn cylinder(
        origin: Point3,
        axis: Vector3,
        radial: Vector3,
        height: f64,
    ) -> NurbsSurface<Vector4> {
        let start = origin + radial;
        let revolved =
            RevolutionSurface::by_revolution(Line(start, start + axis * height), origin, axis);
        NurbsSurface::new(
            revolved
                .try_into_homogeneous_bspline_surface()
                .expect("cylinder net"),
        )
    }

    #[test]
    fn crossing_cylinder_export_is_bounded_and_preserves_guide_direction() {
        let first = cylinder(
            Point3::new(0.0, 0.0, -2.0),
            Vector3::unit_z(),
            Vector3::unit_x() * 0.2,
            4.0,
        );
        let second = cylinder(
            Point3::new(-2.0, 0.0, 0.0),
            Vector3::unit_x(),
            Vector3::unit_y() * 0.5,
            4.0,
        );
        let analytic = |angle: f64, sign: f64| {
            let (sin, cos) = angle.sin_cos();
            let y = 0.2 * sin;
            Point3::new(0.2 * cos, y, sign * (0.25 - y * y).sqrt())
        };
        for (sweep, sign) in [(PI, 1.0), (-PI, 1.0), (PI, -1.0), (-PI, -1.0)] {
            let guide = PolylineCurve(
                (0..=64)
                    .map(|i| analytic(sweep * f64::from(i) / 64.0, sign))
                    .collect(),
            );
            let intersection = IntersectionCurve::new(first.clone(), second.clone(), guide.clone());
            let leader = crossing_cylinder_leader(&intersection, BooleanTolerance::new(1e-5))
                .expect("regular crossing cylinder branch");
            assert_eq!(leader.front(), guide.front());
            assert_eq!(leader.back(), guide.back());
            for i in 0..=1024 {
                let t = f64::from(i) / 1024.0;
                let point = leader.evaluate(t);
                assert!(point.distance(analytic(sweep * t, sign)) <= 2e-9);
            }
        }

        let mut wrong = PolylineCurve(
            (0..=64)
                .map(|i| analytic(PI * f64::from(i) / 64.0, 1.0))
                .collect(),
        );
        wrong[32].z = -wrong[32].z;
        let intersection = IntersectionCurve::new(first.clone(), second, wrong);
        assert!(crossing_cylinder_leader(&intersection, BooleanTolerance::new(1e-5)).is_none());
        let parallel = IntersectionCurve::new(
            first.clone(),
            first,
            PolylineCurve(vec![
                Point3::new(0.2, 0.0, 0.0),
                Point3::new(-0.2, 0.0, 0.0),
            ]),
        );
        assert!(crossing_cylinder_leader(&parallel, BooleanTolerance::new(1e-5)).is_none());
    }
}
