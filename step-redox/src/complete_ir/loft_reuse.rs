//! Proof candidate for translated instances of authored polygonal lofts.
//!
//! Matching world-space points is independent of the source sketch frame and
//! cyclic starting vertex. Correspondence across *every* section is preserved;
//! equal section point sets alone would permit a twisted loft. A candidate is
//! admitted only by the caller's complete original-source boundary proof.

use crate::cad_ir::{CadModel, CadNode, LoftSection, NodeId, RigidTransform};
use crate::math3::{add, distance, sub};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;

const LENGTH: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
type Ring = Vec<[f64; 3]>;
type Section = Vec<Ring>;

struct WorldLoft {
    sections: Vec<Section>,
}

fn world_point(transform: &RigidTransform, point: [f64; 2]) -> [f64; 3] {
    std::array::from_fn(|i| {
        let row = transform.matrix[i];
        row[0].mul_add(point[0], row[1].mul_add(point[1], row[3]))
    })
}

impl WorldLoft {
    fn new(sections: &[LoftSection]) -> Option<Self> {
        let sections = sections
            .iter()
            .map(|section| {
                section
                    .profile
                    .loops
                    .iter()
                    .map(|ring| {
                        let points = ring
                            .polygon_points()?
                            .into_iter()
                            .map(|p| world_point(&section.transform, p))
                            .collect::<Vec<_>>();
                        points
                            .iter()
                            .flatten()
                            .all(|v| v.is_finite())
                            .then_some(points)
                    })
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>()?;
        let counts = sections.first()?.iter().map(Vec::len).collect::<Vec<_>>();
        if sections.len() < 2
            || counts.is_empty()
            || counts.iter().any(|&n| n < 3)
            || sections
                .iter()
                .any(|s| s.iter().map(Vec::len).ne(counts.iter().copied()))
        {
            return None;
        }
        Some(Self { sections })
    }

    fn section(&self, index: usize, reverse: bool) -> &Section {
        &self.sections[if reverse {
            self.sections.len() - 1 - index
        } else {
            index
        }]
    }

    fn aligned_ring(
        &self,
        target: &Self,
        source_ring: usize,
        target_ring: usize,
        correspondence: Correspondence,
    ) -> bool {
        let count = self.sections[0][source_ring].len();
        (0..self.sections.len()).all(|i| {
            let a = &self.sections[i][source_ring];
            let b = &target.section(i, correspondence.reverse_sections)[target_ring];
            a.len() == count
                && b.len() == count
                && a.iter().enumerate().all(|(index, &point)| {
                    let index = if correspondence.reverse_ring {
                        (correspondence.offset + count - index) % count
                    } else {
                        (correspondence.offset + index) % count
                    };
                    distance(add(point, correspondence.translation), b[index]) <= LENGTH
                })
        })
    }

    fn remaining_rings_match(&self, target: &Self, correspondence: Correspondence) -> bool {
        let count = self.sections[0].len();
        let mut used = vec![false; count];
        used[0] = true;
        for source_ring in 1..count {
            let corners = self.sections[0][source_ring].len();
            let matches = (1..count)
                .filter(|&target_ring| {
                    !used[target_ring]
                        && (0..corners).any(|offset| {
                            [false, true].into_iter().any(|reverse_ring| {
                                self.aligned_ring(
                                    target,
                                    source_ring,
                                    target_ring,
                                    Correspondence {
                                        offset,
                                        reverse_ring,
                                        ..correspondence
                                    },
                                )
                            })
                        })
                })
                .collect::<Vec<_>>();
            let [target_ring] = matches.as_slice() else {
                return false;
            };
            used[*target_ring] = true;
        }
        true
    }

    fn translation_to(&self, target: &Self) -> Option<[f64; 3]> {
        if self.sections.len() != target.sections.len()
            || self.sections[0].len() != target.sections[0].len()
            || self.sections[0][0].len() != target.sections[0][0].len()
        {
            return None;
        }
        let anchor = self.sections[0][0][0];
        for reverse_sections in [false, true] {
            let first = target.section(0, reverse_sections);
            for (offset, &p) in first[0].iter().enumerate() {
                let translation = sub(p, anchor);
                for reverse_ring in [false, true] {
                    let correspondence = Correspondence {
                        translation,
                        offset,
                        reverse_ring,
                        reverse_sections,
                    };
                    if self.aligned_ring(target, 0, 0, correspondence)
                        && self.remaining_rings_match(target, correspondence)
                    {
                        return Some(translation);
                    }
                }
            }
        }
        None
    }
}

#[derive(Clone, Copy)]
struct Correspondence {
    translation: [f64; 3],
    offset: usize,
    reverse_ring: bool,
    reverse_sections: bool,
}

/// Share editable loft definitions without changing node IDs or per-instance
/// provenance. This is a geometric proposal, not a certificate; the original
/// source must independently authenticate the result before it can replace
/// the unshared model.
pub(super) fn candidate(model: &CadModel) -> Option<CadModel> {
    let mut exemplars = Vec::<(NodeId, WorldLoft)>::new();
    let mut replacements = Vec::new();
    for (index, node) in model.nodes.iter().enumerate() {
        let CadNode::Loft { sections } = node else {
            continue;
        };
        let Some(witness) = WorldLoft::new(sections) else {
            continue;
        };
        if let Some((child, translation)) = exemplars.iter().find_map(|(id, other)| {
            other
                .translation_to(&witness)
                .map(|translation| (*id, translation))
        }) {
            replacements.push((
                index,
                CadNode::Transform {
                    transform: RigidTransform::translation_mm(translation),
                    child,
                },
            ));
        } else {
            exemplars.push((NodeId(index), witness));
        }
    }
    if replacements.is_empty() {
        return None;
    }
    let mut compact = model.clone();
    for (index, node) in replacements {
        compact.nodes[index] = node;
    }
    Some(compact)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cad_ir::{BooleanOp, Profile2d};
    use anyhow::{Context, Result};

    fn fixture(
        translation: [f64; 3],
        rotate: usize,
        reverse_sections: bool,
        twist: bool,
    ) -> Result<CadNode> {
        let mut sections = Vec::new();
        for (i, points) in [
            [[0., 0.], [4., 0.], [3., 2.], [1., 2.]],
            [[0., 0.], [4., 0.], [2.5, 3.], [1.5, 3.]],
        ]
        .into_iter()
        .enumerate()
        {
            let mut points = points.to_vec();
            points.rotate_left((rotate + usize::from(twist && i == 1)) % 4);
            let profile = Profile2d::polygon(points)?;
            let mut transform = RigidTransform::translation_mm(translation);
            transform.matrix[2][3] += if i == 0 { 0.0 } else { 2.0 };
            sections.push(LoftSection { profile, transform });
        }
        if reverse_sections {
            sections.reverse();
        }
        Ok(CadNode::Loft { sections })
    }

    #[test]
    fn shares_translated_lofts_across_cyclic_order_and_section_direction() -> Result<()> {
        let mut model = CadModel::new();
        let a = model.add_node(fixture([0.; 3], 0, false, false)?);
        let b = model.add_node(fixture([7., -3., 1.], 2, true, false)?);
        let root = model.add_node(CadNode::Boolean {
            op: BooleanOp::Union,
            children: vec![a, b],
        });
        model.add_root(root)?;
        let reused = candidate(&model).context("expected shared loft")?;
        let CadNode::Transform { child, transform } = reused.node(b)? else {
            anyhow::bail!("second loft must instantiate first")
        };
        assert_eq!(*child, a);
        assert!(
            distance(
                [
                    transform.matrix[0][3],
                    transform.matrix[1][3],
                    transform.matrix[2][3]
                ],
                [7., -3., 1.]
            ) <= LENGTH
        );
        assert!(reused.complexity_score(root)? < model.complexity_score(root)?);
        Ok(())
    }

    #[test]
    fn refuses_changed_section_correspondence() -> Result<()> {
        let mut model = CadModel::new();
        model.add_node(fixture([0.; 3], 0, false, false)?);
        model.add_node(fixture([7., -3., 1.], 2, false, true)?);
        assert!(candidate(&model).is_none());
        Ok(())
    }
}
