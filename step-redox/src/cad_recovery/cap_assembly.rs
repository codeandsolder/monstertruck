//! Shared exact source-cap assembly in an optional common rigid frame.
use super::{
    compose_union_models,
    solid::{local_frame_transform, recover_source_cap_extrusion_model},
};
use crate::cad_ir::{CadModel, CadNode, NodeId};
use crate::math3::{add, dot, mul, sub};
use crate::solid_extrusions::RecoveredSourceCapExtrusion;
use anyhow::Result;

#[derive(Clone, Copy, Debug)]
pub(super) struct AssemblyFrame {
    origin: [f64; 3],
    axes: [[f64; 3]; 3],
}

impl AssemblyFrame {
    pub(super) fn centered(
        axes: [[f64; 3]; 3],
        points: impl IntoIterator<Item = [f64; 3]>,
    ) -> Option<Self> {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for point in points {
            for (i, axis) in axes.iter().enumerate() {
                let value = dot(point, *axis);
                if !value.is_finite() {
                    return None;
                }
                min[i] = min[i].min(value);
                max[i] = max[i].max(value);
            }
        }
        if min.iter().any(|x| !x.is_finite()) {
            return None;
        }
        let origin = axes
            .into_iter()
            .zip(min.into_iter().zip(max))
            .fold([0.0; 3], |a, (axis, (lo, hi))| {
                add(a, mul(axis, lo.midpoint(hi)))
            });
        Some(Self { origin, axes })
    }

    fn direction(self, v: [f64; 3]) -> [f64; 3] { self.axes.map(|axis| dot(v, axis)) }

    fn cap(self, cap: &RecoveredSourceCapExtrusion) -> RecoveredSourceCapExtrusion {
        let mut cap = cap.clone();
        cap.origin_mm = self.direction(sub(cap.origin_mm, self.origin));
        cap.x_axis = self.direction(cap.x_axis);
        cap.y_axis = self.direction(cap.y_axis);
        cap.z_axis = self.direction(cap.z_axis);
        cap
    }
}

fn compose_caps(
    caps: &[RecoveredSourceCapExtrusion],
    frame: Option<AssemblyFrame>,
) -> Result<(CadModel, NodeId)> {
    let models = caps
        .iter()
        .map(|cap| {
            frame.map_or_else(
                || recover_source_cap_extrusion_model(cap),
                |frame| recover_source_cap_extrusion_model(&frame.cap(cap)),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let parts = models
        .iter()
        .map(|(model, root)| (model, *root))
        .collect::<Vec<_>>();
    let (mut model, mut root) = compose_union_models(&parts)?;
    if let Some(frame) = frame {
        let [x, y, z] = frame.axes;
        root = model.add_node(CadNode::Transform {
            child: root,
            transform: local_frame_transform(frame.origin, x, y, z),
        });
        model.roots.clear();
        model.add_root(root)?;
    }
    model.validate()?;
    Ok((model, root))
}
/// Replay an exact cap assembly and admit it only after full source proof.
pub(super) fn proven_cap_assembly<K: crate::cad_kernel::CadKernel>(
    caps: &[RecoveredSourceCapExtrusion],
    frame: Option<AssemblyFrame>,
    solid_id: u64,
    entities: &[step_p21::ast::EntityInstance],
    kernel: &K,
) -> Result<Option<(CadModel, NodeId)>> {
    let (model, root) = compose_caps(caps, frame)?;
    let Ok(shape) = kernel.evaluate(&model, root) else {
        return Ok(None);
    };
    if !kernel.summarize(&shape).is_single_solid()
        || !super::rebuilt_kernel_boundary_matches_source(kernel, &shape, solid_id, entities)
            .unwrap_or(false)
    {
        return Ok(None);
    }
    Ok(Some((model, root)))
}
