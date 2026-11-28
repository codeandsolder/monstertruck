//! Exact positive-weight affine circular arcs in homogeneous coordinates.
//! A circle's authored center and basis determine the carrier; nearby
//! topological endpoints must never trigger a new three-point circle fit.

use monstertruck_modeling::{
    AffinePoint, BsplineCurve, InnerProduct, KnotVector, NurbsCurve, Point3, Vector3, Vector4,
};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// Construct a directed circular arc or its nonsingular affine image.
///
/// The supplied radial basis vectors already include scale. Minor arcs use
/// a smooth quadratic parameterization with interior knots; major arcs
/// split at quarter turns to keep all rational weights strictly positive.
pub(super) fn affine_arc(
    origin: Point3,
    x: Vector3,
    y: Vector3,
    sweep: f64,
) -> Option<NurbsCurve<Vector4>> {
    if !sweep.is_finite()
        || sweep == 0.0
        || sweep.abs() > TAU
        || !x.cross(&y).norm_squared().is_finite()
        || x.cross(&y).norm_squared() <= 0.0
    {
        return None;
    }
    if sweep.abs() <= PI {
        let (sin_half, cos_half) = (0.5 * sweep).sin_cos();
        let (sin, cos) = sweep.sin_cos();
        let mut curve = NurbsCurve::new(
            BsplineCurve::new(
                KnotVector::bezier_knot(2),
                vec![
                    (origin.position_vector() + x).push(1.0),
                    (origin.position_vector() * cos_half + x * cos_half + y * sin_half)
                        .extend(cos_half),
                    (origin.position_vector() + x * cos + y * sin).push(1.0),
                ],
            )
            .ok()?,
        );
        // Knot insertion preserves the homogeneous conic exactly, including
        // semicircles whose original middle weight is nearly zero.
        for knot in [0.25, 0.5, 0.75] {
            curve.add_knot(knot);
        }
        return curve
            .control_points()
            .iter()
            .all(|h| h.w > 0.0 && [h.x, h.y, h.z, h.w].into_iter().all(f64::is_finite))
            .then_some(curve);
    }
    let count = (1..=4).find(|n| sweep.abs() <= f64::from(*n) * FRAC_PI_2)?;
    let step = sweep / f64::from(count);
    let weight = (0.5 * step).cos();
    let control = |angle: f64, w: f64| {
        let (sin, cos) = angle.sin_cos();
        let h = origin.position_vector() * w + x * cos + y * sin;
        h.extend(w)
    };
    let mut controls = vec![control(0.0, 1.0)];
    let mut knots = vec![0.0; 3];
    for i in 0..count {
        let start = f64::from(i) * step;
        controls.push(control(step.mul_add(0.5, start), weight));
        controls.push(control(start + step, 1.0));
        let knot = f64::from(i + 1) / f64::from(count);
        knots.extend([knot; 2]);
    }
    knots.push(1.0);
    if controls
        .iter()
        .any(|h| !h.x.is_finite() || !h.y.is_finite() || !h.z.is_finite())
    {
        return None;
    }
    Some(NurbsCurve::new(
        BsplineCurve::new(KnotVector::from(knots), controls).ok()?,
    ))
}
