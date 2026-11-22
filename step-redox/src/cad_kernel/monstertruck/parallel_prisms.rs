//! Boolean algebra of parallel extrusions reduces to planar material algebra
//! on each open interval between cap stations. No volume intersection is
//! necessary; the interval solids meet only at certified separating planes.
use super::{
    BooleanOp, CadModel, CadNode, GeometryPolicy, Matrix4, NodeId, Profile2d, RigidTransform,
    Solid, Surface, Vector3, Vector4, builder, monstertruck_policy, profile, profile_wires,
    rigid_matrix,
};
use anyhow::Result;
use monstertruck_modeling::{Face, SquareMatrix};
use std::collections::HashMap;

struct Prism<'a> {
    profile: &'a Profile2d,
    length: f64,
    transform: Matrix4,
}

enum Instruction {
    Leaf(usize),
    Boolean(BooleanOp, Vec<usize>),
}

#[derive(Default)]
struct Program<'a> {
    prisms: Vec<Prism<'a>>,
    instructions: Vec<Instruction>,
    memo: HashMap<(NodeId, [u64; 16]), usize>,
}

impl<'a> Program<'a> {
    fn emit(&mut self, instruction: Instruction) -> usize {
        let index = self.instructions.len();
        self.instructions.push(instruction);
        index
    }

    fn compile(&mut self, model: &'a CadModel, id: NodeId, matrix: Matrix4) -> Option<usize> {
        let values: &[f64; 16] = matrix.as_ref();
        let key = (id, values.map(f64::to_bits));
        if let Some(&index) = self.memo.get(&key) {
            return Some(index);
        }
        let index = match model.node(id).ok()? {
            CadNode::Transform { transform, child } => {
                self.compile(model, *child, matrix * rigid_matrix(transform))?
            }
            CadNode::Extrude {
                profile,
                vector_mm: [x, y, z],
            } if *x == 0.0 && *y == 0.0 && z.is_finite() && *z != 0.0 => {
                let leaf = self.prisms.len();
                self.prisms.push(Prism {
                    profile,
                    length: *z,
                    transform: matrix,
                });
                self.emit(Instruction::Leaf(leaf))
            }
            CadNode::Boolean { op, children } => {
                let children = children
                    .iter()
                    .map(|&child| self.compile(model, child, matrix))
                    .collect::<Option<Vec<_>>>()?;
                self.emit(Instruction::Boolean(*op, children))
            }
            // In particular, never fuse across a solver-policy override.
            _ => return None,
        };
        self.memo.insert(key, index);
        Some(index)
    }

    fn occupied(&self, active: &[Option<usize>], inside: &[bool]) -> bool {
        let mut values = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            let value = match instruction {
                Instruction::Leaf(leaf) => active[*leaf].is_some_and(|index| inside[index]),
                Instruction::Boolean(op, children) => match op {
                    BooleanOp::Union => children.iter().any(|&id| values[id]),
                    BooleanOp::Intersection => children.iter().all(|&id| values[id]),
                    BooleanOp::Difference => {
                        children.split_first().is_some_and(|(&first, cutters)| {
                            values[first] && !cutters.iter().any(|&id| values[id])
                        })
                    }
                },
            };
            values.push(value);
        }
        values.last().copied().unwrap_or(false)
    }
}

struct PlacedProfile {
    face: Face,
    low: f64,
    high: f64,
}

fn affine(matrix: Matrix4) -> bool {
    let values: &[f64; 16] = matrix.as_ref();
    values.iter().all(|v| v.is_finite())
        && matrix.x.w == 0.0
        && matrix.y.w == 0.0
        && matrix.z.w == 0.0
        && matrix.w.w.to_bits() == 1.0_f64.to_bits()
}

/// The common chart must preserve the complete axial direction and cap
/// planes, not merely project similar vertices. Identical linear transforms
/// use their exact structural equality rather than a noisy inverse product;
/// other orientations are accepted when their relative matrix is exactly
/// block diagonal between the plane and the extrusion axis.
fn relative_chart(frame: Matrix4, inverse: Matrix4, transform: Matrix4) -> Option<Matrix4> {
    if !affine(transform) {
        return None;
    }
    let mut relative = if frame.x == transform.x && frame.y == transform.y && frame.z == transform.z
    {
        Matrix4::identity()
    } else {
        inverse * transform
    };
    relative.w = inverse * (transform.w - frame.w) + Vector4::new(0., 0., 0., 1.);
    if !affine(relative)
        || relative.x.z != 0.0
        || relative.y.z != 0.0
        || relative.z.x != 0.0
        || relative.z.y != 0.0
        || relative.z.z == 0.0
    {
        return None;
    }
    Some(relative)
}

fn place_profiles(prisms: &[Prism<'_>], frame: Matrix4) -> Result<Option<Vec<PlacedProfile>>> {
    if !affine(frame) {
        return Ok(None);
    }
    let Some(inverse) = frame.invert() else {
        return Ok(None);
    };
    let Some(charts) = prisms
        .iter()
        .map(|prism| relative_chart(frame, inverse, prism.transform))
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let mut placed = Vec::with_capacity(prisms.len());
    for (prism, mut chart) in prisms.iter().zip(charts) {
        let start = chart.w.z;
        let end = chart.z.z.mul_add(prism.length, start);
        if !end.is_finite() || end - start == 0.0 {
            return Ok(None);
        }
        chart.w.z = 0.0;
        chart.z.z = 1.0;
        let face = profile::attach_plane_normalized(profile_wires(prism.profile)?)?;
        let mut face = builder::transformed(&face, chart);
        let Surface::Plane(plane) = face.oriented_surface() else {
            return Ok(None);
        };
        if plane.normal().z < 0.0 {
            face.invert();
        }
        placed.push(PlacedProfile {
            face,
            low: start.min(end),
            high: start.max(end),
        });
    }
    Ok(Some(placed))
}

fn section(
    program: &Program<'_>,
    profiles: &[PlacedProfile],
    low: f64,
    high: f64,
    policy: GeometryPolicy,
) -> Result<Option<Solid>> {
    let mut faces = Vec::new();
    let active = profiles
        .iter()
        .map(|profile| {
            if profile.low <= low && high <= profile.high {
                let index = faces.len();
                faces.push(profile.face.clone());
                Some(index)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let regions = match faces.as_slice() {
        [] => Vec::new(),
        [face] => {
            if program.occupied(&active, &[true]) {
                vec![face.clone()]
            } else {
                Vec::new()
            }
        }
        _ => {
            let Some(regions) = monstertruck_solid::try_planar_profile_boolean(
                &faces,
                |inside| program.occupied(&active, inside),
                monstertruck_policy(policy),
            )?
            else {
                return Ok(None);
            };
            regions.into_iter().collect()
        }
    };
    let mut boundaries = Vec::new();
    for face in regions {
        let face =
            builder::transformed(&face, Matrix4::from_translation(Vector3::new(0., 0., low)));
        let body: Solid = builder::extrude(&face, Vector3::new(0., 0., high - low));
        boundaries.extend(body.into_boundaries());
    }
    Ok(Some(Solid::new(boundaries)?))
}

/// Compile the original Boolean DAG once, then evaluate its occupancy in each
/// slab against the same ordered profile leaves. A refused arrangement
/// delegates to the general kernel, without changing the authored CAD IR.
pub(super) fn evaluate(
    model: &CadModel,
    op: BooleanOp,
    children: &[NodeId],
    policy: GeometryPolicy,
) -> Result<Option<Solid>> {
    if matches!(policy, GeometryPolicy::Equivalent) {
        return Ok(None);
    }
    let mut program = Program::default();
    let identity = rigid_matrix(&RigidTransform::identity());
    let Some(children) = children
        .iter()
        .map(|&child| program.compile(model, child, identity))
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    program.emit(Instruction::Boolean(op, children));
    let Some(frame) = program.prisms.first().map(|prism| prism.transform) else {
        return Ok(None);
    };
    let Some(profiles) = place_profiles(&program.prisms, frame)? else {
        return Ok(None);
    };
    let mut levels = profiles
        .iter()
        .flat_map(|profile| [profile.low, profile.high])
        .collect::<Vec<_>>();
    levels.sort_by(f64::total_cmp);
    levels.dedup_by(|a, b| *a - *b == 0.0);
    let mut slabs = Vec::new();
    for limits in levels.windows(2) {
        let Some(slab) = section(&program, &profiles, limits[0], limits[1], policy)? else {
            return Ok(None);
        };
        if !slab.is_empty() {
            slabs.push(slab);
        }
    }
    let joined = match slabs.as_slice() {
        [] => Solid::empty(),
        [slab] => slab.clone(),
        [first, rest @ ..] => {
            let rest = rest.iter().collect::<Vec<_>>();
            let Some(joined) = monstertruck_solid::try_planar_contact_many(
                first,
                &rest,
                monstertruck_solid::PlanarContactOp::Union,
                monstertruck_policy(policy),
            )?
            else {
                return Ok(None);
            };
            joined
        }
    };
    Ok(Some(builder::transformed(&joined, frame)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cad_ir::{Curve2d, ProfileLoop},
        cad_kernel::{
            CadKernel,
            monstertruck::{EvaluatedShape, MonstertruckKernel, extrude_profile_z},
        },
        cad_recovery::rebuilt_kernel_boundary_matches_source,
    };
    use anyhow::Context;

    fn rectangle(low: [f64; 2], high: [f64; 2]) -> Result<Profile2d> {
        Profile2d::polygon(vec![low, [high[0], low[1]], high, [low[0], high[1]]])
    }

    fn disk(center: [f64; 2], radius: f64) -> Profile2d {
        Profile2d {
            loops: vec![ProfileLoop {
                curves: vec![Curve2d::CircleArc {
                    center_mm: center,
                    radius_mm: radius,
                    start_angle_rad: 0.,
                    end_angle_rad: std::f64::consts::TAU,
                }],
            }],
        }
    }

    fn prism(model: &mut CadModel, profile: Profile2d, offset: [f64; 3], length: f64) -> NodeId {
        let child = model.add_node(CadNode::Extrude {
            profile,
            vector_mm: [0., 0., length],
        });
        model.add_node(CadNode::Transform {
            child,
            transform: RigidTransform::translation_mm(offset),
        })
    }

    fn boolean(model: &mut CadModel, op: BooleanOp, children: &[NodeId]) -> NodeId {
        model.add_node(CadNode::Boolean {
            op,
            children: children.to_vec(),
        })
    }

    fn reference_section(profile: &Profile2d, low: f64, high: f64) -> Result<Solid> {
        Ok(builder::transformed(
            &extrude_profile_z(profile, high - low)?,
            Matrix4::from_translation(Vector3::new(0., 0., low)),
        ))
    }

    #[test]
    fn parallel_prism_boolean_preserves_axial_transitions_and_exact_hole() -> Result<()> {
        let mut model = CadModel::new();
        let narrow = rectangle([0., 0.], [6., 4.])?;
        let wide = rectangle([0., 0.], [8., 4.])?;
        let base = prism(&mut model, narrow.clone(), [0., 0., 0.], 8.);
        let addition = prism(
            &mut model,
            rectangle([0., 0.], [4., 4.])?,
            [4., 0., 6.],
            -4.,
        );
        let joined = boolean(&mut model, BooleanOp::Union, &[base, addition]);
        let hole = disk([5., 2.], 0.75);
        let cutter = prism(&mut model, hole.clone(), [0., 0., 3.], 5.);
        let root = boolean(&mut model, BooleanOp::Difference, &[joined, cutter]);
        model.add_root(root)?;
        model.validate()?;
        let policy = GeometryPolicy::Standard {
            tolerance_mm: GeometryPolicy::PROJECT_STANDARD_TOLERANCE_MM,
        };
        let actual = evaluate(&model, BooleanOp::Difference, &[joined, cutter], policy)?
            .context("parallel prism expression was not handled by the planar fast path")?;
        let mut holed = wide.clone();
        holed.loops.extend(hole.loops);
        let expected_slabs = [
            reference_section(&narrow, 0., 2.)?,
            reference_section(&wide, 2., 3.)?,
            reference_section(&holed, 3., 6.)?,
            {
                let mut profile = narrow.clone();
                profile.loops.extend(disk([5., 2.], 0.75).loops);
                reference_section(&profile, 6., 8.)?
            },
        ];
        let expected = monstertruck_solid::try_planar_contact_many(
            &expected_slabs[0],
            &expected_slabs[1..].iter().collect::<Vec<_>>(),
            monstertruck_solid::PlanarContactOp::Union,
            monstertruck_policy(policy),
        )?
        .context("reference contact union failed")?;
        let kernel = MonstertruckKernel::project_standard();
        let reference = kernel.to_step(&EvaluatedShape { solid: expected })?;
        let actual = EvaluatedShape { solid: actual };
        assert!(kernel.summarize(&actual).is_single_solid());
        let exchange =
            monstertruck_io::step::load::ParsedStepSource::parse(reference.as_bytes())?.exchange;
        let section = &exchange.data[0];
        let index = crate::step_graph::build_index(&section.entities);
        let solid_id = section
            .entities
            .iter()
            .find_map(|entry| {
                let record = crate::step_graph::simple_record(entry)?;
                (record.name == "MANIFOLD_SOLID_BREP").then(|| crate::step_graph::entity_id(entry))
            })
            .context("reference STEP omitted solid")?;
        assert!(index.contains_key(&solid_id));
        assert!(rebuilt_kernel_boundary_matches_source(
            &kernel,
            &actual,
            solid_id,
            &section.entities
        )?);
        Ok(())
    }

    #[test]
    fn parallel_prism_intersection_and_empty_difference() -> Result<()> {
        let mut model = CadModel::new();
        let a = prism(&mut model, rectangle([0., 0.], [4., 4.])?, [0., 0., 0.], 4.);
        let b = prism(&mut model, rectangle([0., 0.], [4., 4.])?, [2., 2., 2.], 4.);
        let policy = GeometryPolicy::Standard {
            tolerance_mm: GeometryPolicy::PROJECT_STANDARD_TOLERANCE_MM,
        };
        let intersection = evaluate(&model, BooleanOp::Intersection, &[a, b], policy)?
            .context("parallel intersection fast path refused")?;
        assert_eq!(intersection.boundaries().len(), 1);
        assert_eq!(intersection.boundaries()[0].len(), 6);
        let empty = evaluate(&model, BooleanOp::Difference, &[a, a], policy)?
            .context("coincident difference fast path refused")?;
        assert!(empty.is_empty());
        Ok(())
    }
}
