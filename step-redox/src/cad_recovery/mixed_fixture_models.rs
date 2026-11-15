//! Test-only source-measured analytic features for the 4a6356 engineering
//! fixture. These are *editable CAD sketches and constructive nodes*, not
//! serialized/imported B-reps. They are not a recognition algorithm and must
//! not be used to claim automatic Complete IR admission.
use crate::cad_ir::{
    BooleanOp, CadModel, CadNode, Curve2d, Profile2d, ProfileLoop, RigidTransform,
};
use anyhow::Result;

/// Two original source-measured analytic arc terminal cutters.
pub(crate) fn rounding_caps() -> Result<(CadModel, crate::cad_ir::NodeId, crate::cad_ir::NodeId)> {
    let xc = 67.7545_f64;
    let r = 0.1524_f64;
    let y0 = 5.161_372_032_35_f64;
    let cy = 5.364_572_032_35_f64;
    let cz = 1.679_434_121_45_f64;
    let inner_r = 0.3429_f64;
    let outer_r = 0.6731_f64;
    // Extend the analytical carrier beyond the intersecting feature.
    // Its artificial face boundaries are not part of the target CAD.
    // This also avoids a redundant coplanar interface in the CSG solver.
    let y1 = y0 + 2.0 * r;
    let band_y0 = y0 - 0.03;
    let lower_z = |y: f64| cz + (inner_r.powi(2) - (y - cy).powi(2)).sqrt();
    let upper_z = |y: f64| cz + (outer_r.powi(2) - (y - cy).powi(2)).sqrt();
    let lower_angle = |y: f64| (lower_z(y) - cz).atan2(y - cy);
    let upper_angle = |y: f64| (upper_z(y) - cz).atan2(y - cy);

    let mut model = CadModel::new();
    let band = model.add_node(CadNode::Extrude {
        profile: Profile2d {
            loops: vec![ProfileLoop {
                curves: vec![
                    Curve2d::CircleArc {
                        center_mm: [cy, cz],
                        radius_mm: inner_r,
                        start_angle_rad: lower_angle(band_y0),
                        end_angle_rad: lower_angle(y1),
                    },
                    Curve2d::Line {
                        start_mm: [y1, lower_z(y1)],
                        end_mm: [y1, upper_z(y1)],
                    },
                    Curve2d::CircleArc {
                        center_mm: [cy, cz],
                        radius_mm: outer_r,
                        start_angle_rad: upper_angle(y1),
                        end_angle_rad: upper_angle(band_y0),
                    },
                    Curve2d::Line {
                        start_mm: [band_y0, upper_z(band_y0)],
                        end_mm: [band_y0, lower_z(band_y0)],
                    },
                ],
            }],
        },
        vector_mm: [0.0, 0.0, 2.0 * r + 0.04],
    });
    let band = model.add_node(CadNode::Transform {
        transform: RigidTransform {
            matrix: [
                [0.0, 0.0, 1.0, xc - r - 0.02],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        child: band,
    });
    let top_round = model.add_node(CadNode::Extrude {
        profile: Profile2d {
            loops: vec![ProfileLoop {
                curves: vec![
                    Curve2d::CircleArc {
                        center_mm: [xc, y0],
                        radius_mm: r,
                        start_angle_rad: 0.0,
                        end_angle_rad: std::f64::consts::PI,
                    },
                    Curve2d::Line {
                        start_mm: [xc - r, y0],
                        end_mm: [xc + r, y0],
                    },
                ],
            }],
        },
        vector_mm: [0.0, 0.0, 0.8],
    });
    let top_round = model.add_node(CadNode::Transform {
        transform: RigidTransform {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, cz],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        child: top_round,
    });
    let upper = model.add_node(CadNode::Boolean {
        op: BooleanOp::Intersection,
        children: vec![band, top_round],
    });

    let yb = 6.177_372_032_36_f64;
    let zb = -3.057_665_878_55_f64;
    let bottom = model.add_node(CadNode::Extrude {
        profile: Profile2d {
            loops: vec![ProfileLoop {
                curves: vec![
                    Curve2d::CircleArc {
                        center_mm: [zb, xc],
                        radius_mm: r,
                        start_angle_rad: std::f64::consts::FRAC_PI_2,
                        end_angle_rad: 3.0 * std::f64::consts::FRAC_PI_2,
                    },
                    Curve2d::Line {
                        start_mm: [zb, xc - r],
                        end_mm: [zb, xc + r],
                    },
                ],
            }],
        },
        vector_mm: [0.0, 0.0, 0.3302],
    });
    let lower = model.add_node(CadNode::Transform {
        transform: RigidTransform {
            matrix: [
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, yb],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        child: bottom,
    });
    model.validate()?;
    Ok((model, upper, lower))
}

/// Exact annular arc extrusion intersected with the source-measured taper.
pub(crate) fn curved_lip() -> Result<(CadModel, crate::cad_ir::NodeId)> {
    let x0 = 67.1195;
    let x1 = 68.3895;
    let y0: f64 = 5.351_872_032_35;
    let y1 = 5.859_872_032_35;
    let z0 = 2.022_334_121_45;
    let z1 = 2.352_534_121_45;
    let cy = 5.364_572_032_35;
    let cz = 1.679_434_121_45;
    let lower_r: f64 = 0.3429;
    let upper_r: f64 = 0.6731;
    let lower_at_y0 = cz + (lower_r.powi(2) - (y0 - cy).powi(2)).sqrt();
    let upper_at_y0 = cz + (upper_r.powi(2) - (y0 - cy).powi(2)).sqrt();

    let profile = Profile2d {
        loops: vec![ProfileLoop {
            curves: vec![
                Curve2d::CircleArc {
                    center_mm: [cy, cz],
                    radius_mm: lower_r,
                    start_angle_rad: (lower_at_y0 - cz).atan2(y0 - cy),
                    end_angle_rad: (z0 - cz).atan2(0.0),
                },
                Curve2d::Line {
                    start_mm: [cy, z0],
                    end_mm: [y1, z0],
                },
                Curve2d::Line {
                    start_mm: [y1, z0],
                    end_mm: [y1, z1],
                },
                Curve2d::Line {
                    start_mm: [y1, z1],
                    end_mm: [cy, z1],
                },
                Curve2d::CircleArc {
                    center_mm: [cy, cz],
                    radius_mm: upper_r,
                    start_angle_rad: (z1 - cz).atan2(0.0),
                    end_angle_rad: (upper_at_y0 - cz).atan2(y0 - cy),
                },
                Curve2d::Line {
                    start_mm: [y0, upper_at_y0],
                    end_mm: [y0, lower_at_y0],
                },
            ],
        }],
    };

    let mut model = CadModel::new();
    let arc_body = model.add_node(CadNode::Extrude {
        profile,
        vector_mm: [0.0, 0.0, x1 - x0],
    });
    let arc_world = model.add_node(CadNode::Transform {
        transform: RigidTransform {
            matrix: [
                [0.0, 0.0, 1.0, x0],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        child: arc_body,
    });
    let trapezoid =
        Profile2d::polygon(vec![[x0, y0], [x1, y0], [x1 - 0.254, y1], [x0 + 0.254, y1]])?;
    let prism = model.add_node(CadNode::Extrude {
        profile: trapezoid,
        vector_mm: [0.0, 0.0, z1 - z0 + 0.1],
    });
    let prism_world = model.add_node(CadNode::Transform {
        transform: RigidTransform {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, z0 - 0.05],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        child: prism,
    });
    let root = model.add_node(CadNode::Boolean {
        op: BooleanOp::Intersection,
        children: vec![arc_world, prism_world],
    });
    model.add_root(root)?;
    model.validate()?;
    Ok((model, root))
}
