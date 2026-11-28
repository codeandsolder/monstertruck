//! Join successive regions when their complete ordered section boundaries
//! agree. A station is kept as an editable loft section when the affine
//! trajectories change, rather than becoming two Boolean operands.
use super::super::{extrusion_node, local_frame_transform};
use super::{Frame, IDENTITY_MM, same_profiles};
use crate::{
    cad_ir::{CadModel, CadNode, Curve2d, LoftSection, NodeId, Profile2d, ProfileLoop},
    math2::distance,
    math3::{add, mul},
};
use std::collections::BTreeMap;

#[derive(Clone)]
pub(super) struct Span {
    pub(super) limits: [usize; 2],
    pub(super) start: Profile2d,
    pub(super) end: Profile2d,
}

pub(super) struct Chain {
    profiles: Vec<(usize, Profile2d)>,
}

impl Chain {
    pub(super) fn new(span: Span) -> Self {
        Self {
            profiles: vec![(span.limits[0], span.start), (span.limits[1], span.end)],
        }
    }

    pub(super) fn first_station(&self) -> usize { self.profiles[0].0 }

    fn final_profile(&self) -> &Profile2d { &self.profiles[self.profiles.len() - 1].1 }

    pub(super) fn into_node(
        mut self,
        model: &mut CadModel,
        frame: Frame,
        levels: &[f64],
    ) -> NodeId {
        let lo = levels[self.first_station()];
        let hi = levels[self.profiles[self.profiles.len() - 1].0];
        if self
            .profiles
            .iter()
            .all(|(_, p)| same_profiles(&self.profiles[0].1, p))
        {
            return extrusion_node(model, frame, lo, hi, self.profiles.remove(0).1);
        }
        model.add_node(CadNode::Loft {
            sections: self
                .profiles
                .into_iter()
                .map(|(i, profile)| LoftSection {
                    profile,
                    transform: local_frame_transform(
                        add(frame.origin, mul(frame.axis, levels[i])),
                        frame.x,
                        frame.y,
                        frame.axis,
                    ),
                })
                .collect(),
        })
    }
}

pub(super) fn matching_curve(a: &Curve2d, b: &Curve2d) -> bool {
    match (a, b) {
        (
            Curve2d::Line {
                start_mm: a,
                end_mm: ae,
            },
            Curve2d::Line {
                start_mm: b,
                end_mm: be,
            },
        ) => distance(*a, *b) <= IDENTITY_MM && distance(*ae, *be) <= IDENTITY_MM,
        (Curve2d::CircleArc { .. }, Curve2d::CircleArc { .. }) => a == b,
        _ => false,
    }
}

fn cyclic_offset(a: &ProfileLoop, b: &ProfileLoop) -> Option<usize> {
    let n = a.curves.len();
    if n != b.curves.len() || n == 0 {
        return None;
    }
    let mut matches = (0..n).filter(|&offset| {
        a.curves
            .iter()
            .enumerate()
            .all(|(i, a)| matching_curve(a, &b.curves[(i + offset) % n]))
    });
    let first = matches.next()?;
    matches.next().is_none().then_some(first)
}

/// Match outer ring to outer ring and every hole to one uniquely
/// corresponding hole. Reversal is not allowed: both section boundaries
/// already carry the material orientation.
fn correspondence(a: &Profile2d, b: &Profile2d) -> Option<Vec<(usize, usize)>> {
    if a.loops.len() != b.loops.len() {
        return None;
    }
    let mut used = vec![false; b.loops.len()];
    let mut order = Vec::with_capacity(a.loops.len());
    for (i, ring) in a.loops.iter().enumerate() {
        let mut matches = b.loops.iter().enumerate().filter_map(|(j, other)| {
            if used[j] || (i == 0) != (j == 0) {
                return None;
            }
            Some((j, cyclic_offset(ring, other)?))
        });
        let (j, offset) = matches.next()?;
        if matches.next().is_some() {
            return None;
        }
        used[j] = true;
        order.push((j, offset));
    }
    Some(order)
}

fn align_end(mut end: Profile2d, order: &[(usize, usize)]) -> Profile2d {
    let mut inverse = vec![(0, 0); order.len()];
    for (destination, &(source, rotation)) in order.iter().enumerate() {
        inverse[source] = (destination, rotation);
    }
    let mut rings = end
        .loops
        .into_iter()
        .enumerate()
        .map(|(source, mut ring)| {
            let (destination, rotation) = inverse[source];
            ring.curves.rotate_left(rotation);
            (destination, ring)
        })
        .collect::<Vec<_>>();
    rings.sort_by_key(|(destination, _)| *destination);
    end.loops = rings.into_iter().map(|(_, ring)| ring).collect();
    end
}

/// Only one-to-one section continuations may eliminate an interface.
/// Splitting/joining regions remain independent Boolean operands. The
/// unmodified complete source boundary certificate is still the final
/// admission gate.
pub(super) fn join(spans: Vec<Span>) -> Vec<Chain> {
    let mut starts = BTreeMap::<usize, Vec<Span>>::new();
    for span in spans {
        starts.entry(span.limits[0]).or_default().push(span);
    }
    let mut endings = BTreeMap::<usize, Vec<Option<Chain>>>::new();
    let mut completed = Vec::new();
    for (station, group) in starts {
        let mut previous = endings.remove(&station).unwrap_or_default();
        let candidates = group
            .iter()
            .map(|span| {
                previous
                    .iter()
                    .enumerate()
                    .filter_map(|(i, chain)| {
                        Some((
                            i,
                            correspondence(chain.as_ref()?.final_profile(), &span.start)?,
                        ))
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut counts = vec![0_usize; previous.len()];
        for (i, _) in candidates.iter().flatten() {
            counts[*i] += 1;
        }
        for (span, candidate) in group.into_iter().zip(candidates) {
            let last = span.limits[1];
            let chain = if let [(i, order)] = candidate.as_slice()
                && counts[*i] == 1
                && let Some(mut chain) = previous[*i].take()
            {
                chain.profiles.push((last, align_end(span.end, order)));
                chain
            } else {
                Chain::new(span)
            };
            endings.entry(last).or_default().push(Some(chain));
        }
        completed.extend(previous.into_iter().flatten());
    }
    completed.extend(endings.into_values().flatten().flatten());
    completed.sort_by_key(Chain::first_station);
    completed
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    fn rectangle(w: f64, rotation: usize) -> Result<Profile2d> {
        let mut p = vec![[0., 0.], [w, 0.], [w, 1.], [0., 1.]];
        p.rotate_left(rotation);
        Profile2d::polygon(p)
    }

    #[test]
    fn continuation_preserves_correspondence_and_retains_changed_slope() -> Result<()> {
        let spans = vec![
            Span {
                limits: [0, 1],
                start: rectangle(1., 0)?,
                end: rectangle(2., 0)?,
            },
            Span {
                limits: [1, 2],
                start: rectangle(2., 2)?,
                end: rectangle(4., 2)?,
            },
        ];
        let chains = join(spans);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].profiles.len(), 3);
        assert_eq!(chains[0].final_profile(), &rectangle(4., 0)?);
        Ok(())
    }

    #[test]
    fn ambiguous_continuations_do_not_eliminate_interfaces() -> Result<()> {
        let spans = vec![
            Span {
                limits: [0, 1],
                start: rectangle(1., 0)?,
                end: rectangle(2., 0)?,
            },
            Span {
                limits: [0, 1],
                start: rectangle(3., 0)?,
                end: rectangle(2. + 1.5 * IDENTITY_MM, 0)?,
            },
            Span {
                limits: [1, 2],
                start: rectangle(2., 0)?,
                end: rectangle(3., 0)?,
            },
            Span {
                limits: [1, 2],
                start: rectangle(2. + 0.75 * IDENTITY_MM, 0)?,
                end: rectangle(4., 0)?,
            },
        ];
        assert_eq!(join(spans).len(), 4);
        Ok(())
    }

    #[test]
    fn split_and_gap_keep_independent_regions() -> Result<()> {
        let spans = vec![
            Span {
                limits: [0, 1],
                start: rectangle(1., 0)?,
                end: rectangle(2., 0)?,
            },
            Span {
                limits: [1, 2],
                start: rectangle(3., 0)?,
                end: rectangle(4., 0)?,
            },
            Span {
                limits: [3, 4],
                start: rectangle(4., 0)?,
                end: rectangle(5., 0)?,
            },
        ];
        assert_eq!(join(spans).len(), 3);
        Ok(())
    }
}
