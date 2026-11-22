//! Necessary material-envelope test. Unknown geometry never rejects a proposal.
//! Convex hulls of polynomial and positive rational curves bound the complete
//! sections; affine transforms, extrusion and ruled loft preserve enclosure.
use crate::cad_ir::{
    BooleanOp, CadModel, CadNode, Curve2d, NodeId, Primitive, Profile2d, RigidTransform,
};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug)]
struct Extent {
    min: [f64; 3],
    max: [f64; 3],
}

impl Extent {
    fn new(min: [f64; 3], max: [f64; 3]) -> Option<Self> {
        (0..3)
            .all(|i| min[i].is_finite() && max[i].is_finite() && min[i] <= max[i])
            .then_some(Self { min, max })
    }

    fn point(point: [f64; 3]) -> Option<Self> { Self::new(point, point) }

    fn union(self, other: Self) -> Self {
        Self {
            min: std::array::from_fn(|i| self.min[i].min(other.min[i])),
            max: std::array::from_fn(|i| self.max[i].max(other.max[i])),
        }
    }

    fn transformed(self, transform: &RigidTransform) -> Option<Self> {
        let m = &transform.matrix;
        if m[3].map(f64::to_bits) != [0, 0, 0, 1.0f64.to_bits()]
            || m.iter().flatten().any(|x| !x.is_finite())
        {
            return None;
        }
        let mut min = [0.; 3];
        let mut max = [0.; 3];
        for axis in 0..3 {
            let mut lo = m[axis][3];
            let mut hi = lo;
            for i in 0..3 {
                let a = m[axis][i] * self.min[i];
                let b = m[axis][i] * self.max[i];
                lo = (lo + a.min(b).next_down()).next_down();
                hi = (hi + a.max(b).next_up()).next_up();
            }
            min[axis] = lo;
            max[axis] = hi;
        }
        Self::new(min, max)
    }

    fn contains(self, p: [f64; 3], pad: f64) -> bool {
        (0..3).all(|i| {
            p[i] >= (self.min[i] - pad).next_down() && p[i] <= (self.max[i] + pad).next_up()
        })
    }
}

fn point_hull(points: impl IntoIterator<Item = [f64; 2]>) -> Option<Extent> {
    let mut points = points.into_iter();
    let [x, y] = points.next()?;
    let mut extent = Extent::point([x, y, 0.])?;
    for [x, y] in points {
        extent = extent.union(Extent::point([x, y, 0.])?);
    }
    Some(extent)
}

fn curve_extent(curve: &Curve2d) -> Option<Extent> {
    match curve {
        Curve2d::Line { start_mm, end_mm } => point_hull([*start_mm, *end_mm]),
        Curve2d::CircleArc {
            center_mm: [x, y],
            radius_mm: r,
            ..
        } if *r > 0. => Extent::new(
            [(x - r).next_down(), (y - r).next_down(), 0.],
            [(x + r).next_up(), (y + r).next_up(), 0.],
        ),
        Curve2d::Bezier { control_points_mm } => point_hull(control_points_mm.iter().copied()),
        Curve2d::BSpline {
            control_points_mm,
            weights,
            ..
        } if weights.as_ref().is_none_or(|w| {
            w.len() == control_points_mm.len() && w.iter().all(|x| x.is_finite() && *x > 0.)
        }) =>
        {
            point_hull(control_points_mm.iter().copied())
        }
        _ => None,
    }
}

fn hull(extents: impl IntoIterator<Item = Option<Extent>>) -> Option<Extent> {
    let mut extents = extents.into_iter();
    let mut result = extents.next()??;
    for extent in extents {
        result = result.union(extent?);
    }
    Some(result)
}

fn profile_extent(profile: &Profile2d) -> Option<Extent> {
    hull(
        profile
            .loops
            .iter()
            .flat_map(|l| l.curves.iter().map(curve_extent)),
    )
}

fn compute_extent(
    model: &CadModel,
    id: NodeId,
    cache: &mut HashMap<NodeId, Option<Extent>>,
) -> Option<Extent> {
    if let Some(value) = cache.get(&id) {
        return *value;
    }
    // A defensive unknown entry also terminates a malformed cycle.
    cache.insert(id, None);
    let extent = match model.node(id).ok()? {
        CadNode::Primitive(Primitive::Box { size_mm }) => Extent::new([0.; 3], *size_mm),
        CadNode::Extrude { profile, vector_mm } => {
            let base = profile_extent(profile)?;
            let tip = base.transformed(&RigidTransform::translation_mm(*vector_mm))?;
            Some(base.union(tip))
        }
        CadNode::Loft { sections } => {
            // The IR loft is a ruled interpolation of the placed sections.
            hull(
                sections
                    .iter()
                    .map(|s| profile_extent(&s.profile)?.transformed(&s.transform)),
            )
        }
        CadNode::Transform { child, transform } => {
            compute_extent(model, *child, cache)?.transformed(transform)
        }
        CadNode::Boolean {
            op: BooleanOp::Union,
            children,
        }
        | CadNode::Assembly { children } => {
            hull(children.iter().map(|id| compute_extent(model, *id, cache)))
        }
        CadNode::Boolean {
            op: BooleanOp::Difference | BooleanOp::Intersection,
            children,
        } => compute_extent(model, *children.first()?, cache),
        // Unsupported nodes and explicit solver overrides do not imply any
        // enclosure. Never substitute a sample-based approximation here.
        _ => None,
    };
    cache.insert(id, extent);
    extent
}

/// False proves an original boundary witness lies outside the material
/// enclosure. True only means the proposal remains eligible for full proof.
pub(super) fn covers_witnesses(
    model: &CadModel,
    root: NodeId,
    witnesses: &[[f64; 3]],
    pad: f64,
) -> bool {
    if witnesses.is_empty() || !pad.is_finite() || pad < 0. {
        return true;
    }
    compute_extent(model, root, &mut HashMap::new()).is_none_or(|extent| {
        witnesses
            .iter()
            .all(|&p| p.iter().any(|x| !x.is_finite()) || extent.contains(p, pad))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cad_ir::{LoftSection, ProfileLoop};

    #[test]
    fn transformed_ruled_loft_contains_all_rational_control_extrema() {
        let profile = Profile2d {
            loops: vec![ProfileLoop {
                curves: vec![Curve2d::BSpline {
                    degree: 2,
                    control_points_mm: vec![[0., 0.], [10., 3.], [2., 0.]],
                    knots: vec![0., 0., 0., 1., 1., 1.],
                    weights: Some(vec![1., 2., 1.]),
                }],
            }],
        };
        let mut m = CadModel::new();
        let root = m.add_node(CadNode::Loft {
            sections: vec![
                LoftSection {
                    profile: profile.clone(),
                    transform: RigidTransform::identity(),
                },
                LoftSection {
                    profile,
                    transform: RigidTransform::translation_mm([-3., 0., 5.]),
                },
            ],
        });
        assert!(covers_witnesses(
            &m,
            root,
            &[[10., 3., 0.], [-3., 0., 5.]],
            0.
        ));
        assert!(!covers_witnesses(&m, root, &[[10.1, 3., 0.]], 0.));
    }

    #[test]
    fn unsupported_geometry_does_not_prune_union() {
        let mut m = CadModel::new();
        let b = m.add_node(CadNode::Primitive(Primitive::Box { size_mm: [1.; 3] }));
        let s = m.add_node(CadNode::Primitive(Primitive::Sphere { radius_mm: 10. }));
        let union = m.add_node(CadNode::Boolean {
            op: BooleanOp::Union,
            children: vec![b, s],
        });
        assert!(covers_witnesses(&m, union, &[[8., 0., 0.]], 0.));
        let difference = m.add_node(CadNode::Boolean {
            op: BooleanOp::Difference,
            children: vec![b, s],
        });
        assert!(!covers_witnesses(&m, difference, &[[8., 0., 0.]], 0.));
        assert!(covers_witnesses(&m, b, &[[1.0001, 0., 0.]], 0.001));
    }
}
