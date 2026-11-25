//! Necessary material-envelope test. Unknown geometry never rejects a proposal.
//! Convex hulls of polynomial and positive rational curves bound the complete
//! sections; affine transforms, extrusion and ruled loft preserve enclosure.
use crate::cad_ir::{
    BooleanOp, CadModel, CadNode, Curve2d, NodeId, Primitive, Profile2d, RigidTransform,
};
use std::collections::{HashMap, HashSet};

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

/// A necessary occupancy expression, not a geometric classification.
/// Union preserves the separate operand envelopes; intersection requires
/// every known envelope; subtraction can only occupy its positive operand.
enum BoundTest {
    Unknown,
    Extent(Extent),
    Any(Vec<usize>),
    All(Vec<usize>),
}

#[derive(Default)]
struct WitnessProgram {
    tests: Vec<BoundTest>,
    memo: HashMap<(NodeId, Vec<NodeId>), usize>,
    visiting: HashSet<NodeId>,
}

impl WitnessProgram {
    fn push(&mut self, test: BoundTest) -> usize {
        let id = self.tests.len();
        self.tests.push(test);
        id
    }

    fn compile(&mut self, model: &CadModel, id: NodeId, transforms: &mut Vec<NodeId>) -> usize {
        let key = (id, transforms.clone());
        if let Some(&test) = self.memo.get(&key) {
            return test;
        }
        if !self.visiting.insert(id) {
            // Invalid cyclic model: no geometric rejection is justified.
            return self.push(BoundTest::Unknown);
        }
        let index = match model.node(id) {
            Ok(CadNode::Transform { child, .. }) => {
                transforms.push(id);
                let child = self.compile(model, *child, transforms);
                transforms.pop();
                child
            }
            Ok(CadNode::Boolean {
                op: BooleanOp::Difference,
                children,
            }) => match children.first() {
                Some(&child) => self.compile(model, child, transforms),
                None => self.push(BoundTest::Unknown),
            },
            Ok(
                CadNode::Boolean {
                    op: BooleanOp::Union,
                    children,
                }
                | CadNode::Assembly { children },
            ) => {
                let operands = children
                    .iter()
                    .map(|&child| self.compile(model, child, transforms))
                    .collect();
                self.push(BoundTest::Any(operands))
            }
            Ok(CadNode::Boolean {
                op: BooleanOp::Intersection,
                children,
            }) => {
                let operands = children
                    .iter()
                    .map(|&child| self.compile(model, child, transforms))
                    .collect();
                self.push(BoundTest::All(operands))
            }
            Ok(node) => {
                let extent = leaf_extent(node).and_then(|mut extent| {
                    // Apply each original affine map with outward rounding;
                    // matrix composition would lose that enclosure guarantee.
                    for &transform in transforms.iter().rev() {
                        let CadNode::Transform { transform, .. } = model.node(transform).ok()?
                        else {
                            return None;
                        };
                        extent = extent.transformed(transform)?;
                    }
                    Some(extent)
                });
                self.push(extent.map_or(BoundTest::Unknown, BoundTest::Extent))
            }
            Err(_) => self.push(BoundTest::Unknown),
        };
        self.visiting.remove(&id);
        self.memo.insert(key, index);
        index
    }

    fn covers(&self, root: usize, point: [f64; 3], pad: f64, values: &mut Vec<bool>) -> bool {
        values.clear();
        for test in &self.tests {
            values.push(match test {
                BoundTest::Unknown => true,
                BoundTest::Extent(extent) => extent.contains(point, pad),
                BoundTest::Any(operands) => {
                    operands.is_empty() || operands.iter().any(|&id| values[id])
                }
                BoundTest::All(operands) => operands.iter().all(|&id| values[id]),
            });
        }
        values[root]
    }
}

fn leaf_extent(node: &CadNode) -> Option<Extent> {
    match node {
        CadNode::Primitive(Primitive::Box { size_mm }) => Extent::new([0.; 3], *size_mm),
        CadNode::Extrude { profile, vector_mm } => {
            let base = profile_extent(profile)?;
            let tip = base.transformed(&RigidTransform::translation_mm(*vector_mm))?;
            Some(base.union(tip))
        }
        CadNode::Loft { sections } => {
            // Ruled interpolation remains within the hull of placed sections.
            hull(
                sections
                    .iter()
                    .map(|s| profile_extent(&s.profile)?.transformed(&s.transform)),
            )
        }
        // Unknown geometry and tolerance overrides never reject a proposal.
        _ => None,
    }
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
    let mut program = WitnessProgram::default();
    let root = program.compile(model, root, &mut Vec::new());
    let mut values = Vec::with_capacity(program.tests.len());
    witnesses.iter().all(|&point| {
        point.iter().any(|x| !x.is_finite()) || program.covers(root, point, pad, &mut values)
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
    fn disjoint_union_preserves_gaps_through_outer_transforms() {
        let mut m = CadModel::new();
        let b = m.add_node(CadNode::Primitive(Primitive::Box { size_mm: [1.; 3] }));
        let shifted = m.add_node(CadNode::Transform {
            child: b,
            transform: RigidTransform::translation_mm([4., 0., 0.]),
        });
        let union = m.add_node(CadNode::Boolean {
            op: BooleanOp::Union,
            children: vec![b, shifted],
        });
        let moved = m.add_node(CadNode::Transform {
            child: union,
            transform: RigidTransform {
                matrix: [
                    [0., -1., 0., 2.],
                    [1., 0., 0., 3.],
                    [0., 0., 1., 4.],
                    [0., 0., 0., 1.],
                ],
            },
        });
        assert!(covers_witnesses(
            &m,
            moved,
            &[[1.5, 3.5, 4.5], [1.5, 7.5, 4.5]],
            0.
        ));
        assert!(!covers_witnesses(&m, moved, &[[1.5, 5.5, 4.5]], 0.));
        assert!(covers_witnesses(&m, moved, &[[1.5, 4.0001, 4.5]], 0.001));
    }

    #[test]
    fn intersections_require_known_operands_and_differences_only_positive() {
        let mut m = CadModel::new();
        let a = m.add_node(CadNode::Primitive(Primitive::Box { size_mm: [4.; 3] }));
        let b = m.add_node(CadNode::Primitive(Primitive::Box { size_mm: [1.; 3] }));
        let unknown = m.add_node(CadNode::Primitive(Primitive::Sphere { radius_mm: 10. }));
        let intersect = m.add_node(CadNode::Boolean {
            op: BooleanOp::Intersection,
            children: vec![unknown, a, b],
        });
        assert!(covers_witnesses(&m, intersect, &[[0.5; 3]], 0.));
        assert!(!covers_witnesses(&m, intersect, &[[2.; 3]], 0.));
        let diff = m.add_node(CadNode::Boolean {
            op: BooleanOp::Difference,
            children: vec![a, b],
        });
        assert!(covers_witnesses(&m, diff, &[[0.5; 3]], 0.));
        assert!(covers_witnesses(&m, diff, &[[2.; 3]], 0.));
        assert!(!covers_witnesses(&m, diff, &[[5.; 3]], 0.));
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
