//! Reconstruct a constant extrusion host with material changes at witnessed
//! transverse stations.  Unlike a slab union, this preserves continuous
//! analytic side faces across station boundaries.
use super::{DIRECTION, Frame, LENGTH, extrusion_node, section_profiles};
use crate::brep::SurfaceSupport;
use crate::cad_ir::{BooleanOp, CadModel, CadNode, NodeId, Profile2d};
use crate::math2::distance;
use crate::math3::{cross, dot, norm, normalize};
use crate::profile_curves::RecoveredProfileCurve;
use crate::reduced_brep::ReducedSolid;
use crate::solid_extrusions::source_cap_loops_in_frame;
use anyhow::Result;
use std::collections::HashMap;
use step_p21::ast::EntityInstance;

#[derive(Clone)]
struct Event {
    face_id: u64,
    station_index: usize,
    removes: bool,
    loops: Vec<Vec<RecoveredProfileCurve>>,
}

fn same_curve(first: &RecoveredProfileCurve, second: &RecoveredProfileCurve) -> bool {
    match (first, second) {
        (
            RecoveredProfileCurve::Line {
                start_mm: first_start,
                end_mm: first_end,
                ..
            },
            RecoveredProfileCurve::Line {
                start_mm: second_start,
                end_mm: second_end,
                ..
            },
        ) => {
            distance(*first_start, *second_start) <= LENGTH
                && distance(*first_end, *second_end) <= LENGTH
        }
        (
            RecoveredProfileCurve::CircleArc {
                center_mm: first_center,
                radius_mm: first_radius,
                start_angle_rad: first_start,
                end_angle_rad: first_end,
                ..
            },
            RecoveredProfileCurve::CircleArc {
                center_mm: second_center,
                radius_mm: second_radius,
                start_angle_rad: second_start,
                end_angle_rad: second_end,
                ..
            },
        ) => {
            distance(*first_center, *second_center) <= LENGTH
                && (first_radius - second_radius).abs() <= LENGTH
                && ((first_end - first_start) - (second_end - second_start)).abs()
                    * first_radius.max(*second_radius)
                    <= LENGTH
                && first
                    .start_point()
                    .zip(second.start_point())
                    .is_some_and(|(left, right)| distance(left, right) <= LENGTH)
                && first
                    .end_point()
                    .zip(second.end_point())
                    .is_some_and(|(left, right)| distance(left, right) <= LENGTH)
        }
        _ => false,
    }
}

fn same_ring(first: &[RecoveredProfileCurve], second: &[RecoveredProfileCurve]) -> bool {
    first.len() == second.len()
        && (0..first.len()).any(|offset| {
            first
                .iter()
                .zip(second.iter().cycle().skip(offset))
                .all(|(a, b)| same_curve(a, b))
        })
}

fn same_region(first: &Event, second: &Event) -> bool {
    if first.loops.len() != second.loops.len() {
        return false;
    }
    let mut used = vec![false; second.loops.len()];
    for ring in &first.loops {
        let matches = second
            .loops
            .iter()
            .enumerate()
            .filter(|(i, other)| !used[*i] && same_ring(ring, other))
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        let [index] = matches.as_slice() else {
            return false;
        };
        used[*index] = true;
    }
    true
}

fn events(
    solid: &ReducedSolid,
    frame: Frame,
    levels: &[f64],
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<Vec<Event>> {
    let mut events = Vec::new();
    for face in solid.faces() {
        let SurfaceSupport::Plane(plane) = face.surface else {
            continue;
        };
        let normal = normalize(plane.normal, DIRECTION)?;
        if norm(cross(normal, frame.axis)) > DIRECTION {
            continue;
        }
        let at = frame.axial(plane.origin_mm);
        let station_index = levels
            .iter()
            .position(|level| (at - level).abs() <= LENGTH)?;
        if station_index == 0 || station_index + 1 == levels.len() {
            continue;
        }
        let removes = (dot(normal, frame.axis) > 0.0) == face.surface_forward;
        let loops = source_cap_loops_in_frame(
            face.id,
            frame.origin,
            [frame.x, frame.y, frame.axis],
            entities,
            index,
        )?;
        events.push(Event {
            face_id: face.id,
            station_index,
            removes,
            loops,
        });
    }
    events.sort_by_key(|event| (event.station_index, event.face_id));
    (!events.is_empty()).then_some(events)
}

/// Every interior transverse face is a signed occupancy change. Extrude the
/// initial material section for the entire axial span, then replay those
/// changes. Opposite matching faces at adjacent stations form a single
/// bounded change rather than a remove-and-restore pair of suffix operations.
/// Final source-boundary equality is independently checked by the caller.
pub(super) fn model(
    solid: &ReducedSolid,
    frame: Frame,
    levels: &[f64],
    initial: &[Profile2d],
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Result<Option<(CadModel, NodeId)>> {
    let Some(events) = events(solid, frame, levels, entities, index) else {
        return Ok(None);
    };
    let mut model = CadModel::new();
    let mut roots = Vec::new();
    let [first, .., last] = levels else {
        return Ok(None);
    };
    let (first, last) = (*first, *last);
    for profile in initial {
        roots.push(extrusion_node(
            &mut model,
            frame,
            first,
            last,
            profile.clone(),
        ));
    }
    if roots.is_empty() {
        return Ok(None);
    }
    let mut root = super::add_bounded_union_tree(&mut model, roots);
    let mut consumed = vec![false; events.len()];
    for (i, event) in events.iter().enumerate() {
        if consumed[i] {
            continue;
        }
        let matches = events
            .iter()
            .enumerate()
            .skip(i + 1)
            .filter(|(j, other)| {
                !consumed[*j]
                    && other.station_index == event.station_index + 1
                    && event.removes != other.removes
                    && same_region(event, other)
            })
            .map(|(j, _)| j)
            .collect::<Vec<_>>();
        let end = if let [j] = matches.as_slice() {
            // The inverse event at the next station is paired only if no
            // competing event claims that same region. Same-station faces
            // are disjoint material patches of the original boundary.
            let inverse = &events[*j];
            if events[..*j]
                .iter()
                .enumerate()
                .filter(|(k, other)| {
                    !consumed[*k]
                        && other.station_index == event.station_index
                        && other.removes != inverse.removes
                        && same_region(other, inverse)
                })
                .count()
                == 1
            {
                consumed[*j] = true;
                levels[inverse.station_index]
            } else {
                last
            }
        } else {
            last
        };
        let Some(profiles) = section_profiles(&event.loops) else {
            return Ok(None);
        };
        for profile in profiles {
            let change =
                extrusion_node(&mut model, frame, levels[event.station_index], end, profile);
            root = model.add_node(CadNode::Boolean {
                op: if event.removes {
                    BooleanOp::Difference
                } else {
                    BooleanOp::Union
                },
                children: vec![root, change],
            });
        }
    }
    model.add_root(root)?;
    model.validate()?;
    Ok(Some((model, root)))
}
