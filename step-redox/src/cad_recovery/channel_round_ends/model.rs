//! Construct exact semicircular prisms and their coaxial-band intersections.

use super::{EndCarriers, LENGTH, RoundEndProposal, RoundFrame};
use crate::{
    cad_ir::{BooleanOp, CadModel, CadNode, Curve2d, NodeId, ProofStatus, Provenance},
    cad_recovery::solid::extrude_profile_in_frame as extrude,
    math3::{add, mul},
};
use anyhow::{Result, ensure};
use std::f64::consts::PI;

fn semicircle(model: &mut CadModel, frame: RoundFrame, lower: f64, upper: f64) -> NodeId {
    let r = frame.radius;
    let (start_angle, end_angle, start, end) = if frame.positive_side {
        (0.0, PI, [r, 0.0], [-r, 0.0])
    } else {
        (PI, 2.0 * PI, [-r, 0.0], [r, 0.0])
    };
    extrude(
        model,
        vec![
            Curve2d::CircleArc {
                center_mm: [0.0, 0.0],
                radius_mm: r,
                start_angle_rad: start_angle,
                end_angle_rad: end_angle,
            },
            Curve2d::Line {
                start_mm: end,
                end_mm: start,
            },
        ],
        upper - lower,
        add(frame.origin, mul(frame.axis, lower)),
        [frame.width_axis, frame.side_axis, frame.axis],
    )
}

struct CoaxialBand {
    center: [f64; 2],
    inner: f64,
    outer: f64,
    positive_root: bool,
}

impl CoaxialBand {
    fn at(&self, side: f64, radius: f64) -> [f64; 2] {
        let [cy, cz] = self.center;
        let delta = side - cy;
        let axial = ((radius - delta.abs()) * (radius + delta.abs())).sqrt();
        [side, cz + if self.positive_root { axial } else { -axial }]
    }

    fn arc(&self, radius: f64, start: [f64; 2], end: [f64; 2]) -> Curve2d {
        let angle = |p: [f64; 2]| (p[1] - self.center[1]).atan2(p[0] - self.center[0]);
        Curve2d::CircleArc {
            center_mm: self.center,
            radius_mm: radius,
            start_angle_rad: angle(start),
            end_angle_rad: angle(end),
        }
    }

    fn profile(&self, lower: f64, upper: f64) -> Vec<Curve2d> {
        // A negative branch reverses the radial-plane winding.
        let (lo, hi) = if self.positive_root {
            (lower, upper)
        } else {
            (upper, lower)
        };
        let a = self.at(lo, self.inner);
        let b = self.at(hi, self.inner);
        let c = self.at(hi, self.outer);
        let d = self.at(lo, self.outer);
        vec![
            self.arc(self.inner, a, b),
            Curve2d::Line {
                start_mm: b,
                end_mm: c,
            },
            self.arc(self.outer, c, d),
            Curve2d::Line {
                start_mm: d,
                end_mm: a,
            },
        ]
    }

    fn intersection(&self, model: &mut CadModel, frame: RoundFrame) -> Result<NodeId> {
        let signed_r = if frame.positive_side {
            frame.radius
        } else {
            -frame.radius
        };
        let (lo, hi) = (signed_r.min(0.0), signed_r.max(0.0));
        let gap = self.inner - (lo - self.center[0]).abs().max((hi - self.center[0]).abs());
        ensure!(
            gap > LENGTH,
            "semicircle reaches a limiting-cylinder tangency"
        );
        // The extra planes are outside the semicircle's full footprint, and
        // their offsets are bounded by the separated inner-root domain.
        // They are construction boundaries, not fitted source dimensions.
        let clearance = 0.25 * gap.min(frame.radius);
        let extent = frame.radius + clearance;
        let band = extrude(
            model,
            self.profile(lo - clearance, hi + clearance),
            2.0 * extent,
            add(frame.origin, mul(frame.width_axis, -extent)),
            [frame.side_axis, frame.axis, frame.width_axis],
        );
        let height = self.outer + clearance;
        let center = self.center[1];
        let (lower, upper) = if self.positive_root {
            (center, center + height)
        } else {
            (center - height, center)
        };
        let round = semicircle(model, frame, lower, upper);
        Ok(model.add_node(CadNode::Boolean {
            op: BooleanOp::Intersection,
            children: vec![band, round],
        }))
    }
}

impl RoundEndProposal {
    pub(in crate::cad_recovery) fn model(&self) -> Result<(CadModel, NodeId)> {
        let mut model = CadModel::new();
        let root = match self.ends {
            EndCarriers::Planar { lower, upper } => {
                semicircle(&mut model, self.frame, lower, upper)
            }
            EndCarriers::Coaxial {
                center_side,
                center_axial,
                inner_radius,
                outer_radius,
                positive_root,
            } => CoaxialBand {
                center: [center_side, center_axial],
                inner: inner_radius,
                outer: outer_radius,
                positive_root,
            }
            .intersection(&mut model, self.frame)?,
        };
        model.set_provenance(
            root,
            Provenance {
                source_entity_ids: self.source_entity_ids.clone(),
                proof: ProofStatus::HeuristicCandidate,
                max_residual_mm: None,
            },
        )?;
        model.add_root(root)?;
        model.validate()?;
        Ok((model, root))
    }
}
