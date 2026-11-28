use super::*;

#[test]
fn degree_elevated_rational_circle_is_certified_on_the_entire_domain() {
    // Degree elevation in homogeneous coordinates of the unit quarter-circle.
    let w = 0.5_f64.sqrt();
    let q = (1.0 + 2.0 * w) / 3.0;
    let inner = 2.0 * w / (1.0 + 2.0 * w);
    let mut profile = BSplineSupport {
        degree: 3,
        control_points_mm: vec![
            [1.0, 0.0, 0.0],
            [1.0, inner, 0.0],
            [inner, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        weights: Some(vec![1.0, q, q, 1.0]),
    };
    let certify = |spline: &BSplineSupport| {
        certify_rational_bspline_circle(spline, [0.0; 3], [0.0, 0.0, 1.0], 1.0, 1.0e-7)
    };
    assert_eq!(certify(&profile), Some(true));
    let surface = SplineExtrusionSupport {
        profile: profile.clone(),
        extrusion_mm: [0.0, 0.0, 2.0],
        max_residual_mm: 0.0,
    };
    assert!(matches!(
        spline_extrusion_cylinder_key(&surface, true, &mut CellWitnesses::default()),
        Some(OrientedSupportKey::Cylinder {
            radial_outward: true,
            ..
        })
    ));
    profile.control_points_mm[1][2] = 0.001;
    assert_eq!(certify(&profile), Some(false));
    profile.control_points_mm[1][2] = 0.0;
    profile.control_points_mm[1][0] += 0.001;
    assert_eq!(certify(&profile), Some(false));
}
