use super::*;
use monstertruck_io::step::load::step_geometry::{KnotVector, NurbsSurface};

#[test]
fn native_classifier_recovers_both_axes_and_nonunit_affine_spline_charts() -> Result<()> {
    let profile = KnotVector::bezier_knot(2);
    let sweep = KnotVector::try_from(vec![2.0, 2.0, 2.0, 3.0, 5.0, 6.0, 6.0, 6.0])?;
    let coordinates = [0.0, 0.125, 0.5, 0.875, 1.0];
    let points = [(0.0, 0.0, 1.0), (1.0, 1.0, 0.75), (2.0, 0.0, 1.0)];
    for swapped in [false, true] {
        let pole = |i: usize, j: usize| {
            let (x, y, w) = points[i];
            Vector4::new(x * w, y * w, 4.0 * coordinates[j] * w, w)
        };
        let (knots, control) = if swapped {
            (
                (sweep.clone(), profile.clone()),
                (0..5)
                    .map(|j| (0..3).map(|i| pole(i, j)).collect())
                    .collect(),
            )
        } else {
            (
                (profile.clone(), sweep.clone()),
                (0..3)
                    .map(|i| (0..5).map(|j| pole(i, j)).collect())
                    .collect(),
            )
        };
        let net = BsplineSurface::new(knots, control)?;
        let SurfaceSupport::SplineExtrusion(evidence) =
            native_surface_support(&Surface::Nurbs(NurbsSurface::new(net)), 1)
        else {
            anyhow::bail!("native affine spline wall was not recognized");
        };
        let sign = if swapped { -1.0 } else { 1.0 };
        assert_eq!(evidence.extrusion_mm, [0.0, 0.0, sign * 4.0]);
        assert_eq!(evidence.profile.weights, Some(vec![1.0, 0.75, 1.0]));
        assert_eq!(
            evidence.profile.control_points_mm[0],
            [0.0, 0.0, if swapped { 4.0 } else { 0.0 }]
        );
        assert!(evidence.max_residual_mm < 1.0e-12);
    }
    Ok(())
}
