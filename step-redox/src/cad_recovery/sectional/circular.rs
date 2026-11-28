//! Exact interior sections of a cylinder with rectilinear axial/angular trims.
//!
//! On the periodic (angle, axial station) chart, generator edges are vertical
//! and circular edges horizontal. Occupancy at an interior angle is the odd
//! parity of transverse arc crossings below the slice. This distinguishes
//! major/minor sectors, holes, disconnected sectors, and full periodic faces
//! without assuming that every cylinder has a four-edge rectangular trim.
use super::{DIRECTION, Frame, IDENTITY_MM, SliceGraph, valid_line};
use crate::brep::{CurveSupport, OrientedEdgeUse, SurfaceSupport};
use crate::math2::{distance as distance2, sub as sub2};
use crate::math3::{cross, norm, normalize, point_to_unit_line_distance};
use crate::profile_curves::RecoveredProfileCurve;
use crate::reduced_brep::ReducedFace;
use crate::solid_extrusions::{angle_on_arc, recovered_profile_curve};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

#[derive(Debug, Clone, Copy)]
pub(super) struct TransverseArc {
    center: [f64; 2],
    radius: f64,
    start: f64,
    end: f64,
    station: f64,
}

/// Certify the complete source circle locus, including bounded edge sense.
/// A geometrically coincident pair of *distinct* vertices is not a full circle.
pub(super) fn transverse_arc(edge: &OrientedEdgeUse, frame: Frame) -> Option<TransverseArc> {
    let CurveSupport::Circle(circle) = edge.support else {
        return None;
    };
    let normal = normalize(circle.normal, DIRECTION)?;
    if norm(cross(normal, frame.axis)) > DIRECTION
        || !circle.radius_mm.is_finite()
        || circle.radius_mm <= IDENTITY_MM
        || !circle.center_mm.iter().all(|v| v.is_finite())
    {
        return None;
    }
    let center = frame.project(circle.center_mm);
    let station = frame.axial(circle.center_mm);
    for point in [edge.start_mm, edge.end_mm] {
        if (frame.axial(point) - station).abs() > frame.tolerance_mm
            || (distance2(frame.project(point), center) - circle.radius_mm).abs()
                > frame.tolerance_mm
            || !point.iter().all(|v| v.is_finite())
        {
            return None;
        }
    }
    let closed = edge.start_vertex == edge.end_vertex;
    if closed
        != (distance2(frame.project(edge.start_mm), frame.project(edge.end_mm)) <= IDENTITY_MM)
    {
        return None;
    }
    let RecoveredProfileCurve::CircleArc {
        start_angle_rad: start,
        end_angle_rad: end,
        ..
    } = recovered_profile_curve(
        edge,
        frame.origin,
        frame.x,
        frame.y,
        frame.axis,
        frame.tolerance_mm,
    )?
    else {
        return None;
    };
    (start.is_finite() && end.is_finite() && (end - start).abs() * circle.radius_mm > IDENTITY_MM)
        .then_some(TransverseArc {
            center,
            radius: circle.radius_mm,
            start,
            end,
            station,
        })
}

#[derive(Clone, Copy)]
struct Generator {
    id: u64,
    angle: f64,
    low: f64,
    high: f64,
    multiplicity: usize,
}

pub(super) struct CylinderTrim {
    center: [f64; 2],
    radius: f64,
    outward: bool,
    generators: Vec<Generator>,
    arcs: Vec<TransverseArc>,
}

impl CylinderTrim {
    pub(super) fn new(face: &ReducedFace, frame: Frame) -> Option<Self> {
        let SurfaceSupport::Cylinder(cylinder) = face.surface else {
            return None;
        };
        let axis = normalize(cylinder.axis, DIRECTION)?;
        if norm(cross(axis, frame.axis)) > DIRECTION
            || !cylinder.radius_mm.is_finite()
            || cylinder.radius_mm <= IDENTITY_MM
            || !cylinder.axis_origin_mm.iter().all(|v| v.is_finite())
        {
            return None;
        }
        let center = frame.project(cylinder.axis_origin_mm);
        let radius = cylinder.radius_mm;
        let mut generators = BTreeMap::<u64, Generator>::new();
        let mut arcs = Vec::new();
        for edge in face.loop_edges() {
            match edge.support {
                CurveSupport::Circle(circle) => {
                    let arc = transverse_arc(edge, frame)?;
                    if distance2(arc.center, center) > frame.tolerance_mm
                        || (arc.radius - radius).abs() > frame.tolerance_mm
                        || point_to_unit_line_distance(
                            circle.center_mm,
                            cylinder.axis_origin_mm,
                            axis,
                        ) > frame.tolerance_mm
                    {
                        return None;
                    }
                    arcs.push(arc);
                }
                CurveSupport::Line(line) => {
                    if !valid_line(edge, frame.tolerance_mm)
                        || norm(cross(normalize(line.direction, DIRECTION)?, frame.axis))
                            > DIRECTION
                    {
                        return None;
                    }
                    let a = frame.project(edge.start_mm);
                    let b = frame.project(edge.end_mm);
                    if distance2(a, b) > frame.tolerance_mm
                        || (distance2(a, center) - radius).abs() > frame.tolerance_mm
                        || (distance2(b, center) - radius).abs() > frame.tolerance_mm
                    {
                        return None;
                    }
                    let start = frame.axial(edge.start_mm);
                    let end = frame.axial(edge.end_mm);
                    let low = start.min(end);
                    let high = start.max(end);
                    if high - low <= IDENTITY_MM {
                        return None;
                    }
                    let delta = sub2(a, center);
                    let angle = delta[1].atan2(delta[0]).rem_euclid(TAU);
                    if let Some(previous) = generators.get_mut(&edge.edge_id) {
                        if previous.multiplicity != 1
                            || (previous.low - low).abs() > IDENTITY_MM
                            || (previous.high - high).abs() > IDENTITY_MM
                            || distance2(
                                point(center, radius, previous.angle),
                                point(center, radius, angle),
                            ) > IDENTITY_MM
                        {
                            return None;
                        }
                        previous.multiplicity += 1;
                    } else {
                        generators.insert(
                            edge.edge_id,
                            Generator {
                                id: edge.edge_id,
                                angle,
                                low,
                                high,
                                multiplicity: 1,
                            },
                        );
                    }
                }
                _ => return None,
            }
        }
        if arcs.is_empty() {
            return None;
        }
        Some(Self {
            center,
            radius,
            outward: face.surface_forward,
            generators: generators.into_values().collect(),
            arcs,
        })
    }

    /// Check each angular cell, including source endpoints at stations other
    /// than the slice. No internal occupancy transition can be hidden within
    /// an emitted arc between two active generators.
    fn occupied(&self, start: f64, end: f64, level: f64) -> Option<bool> {
        let mut cuts = vec![start, end];
        for arc in &self.arcs {
            for angle in [arc.start, arc.end] {
                let unwrapped = start + (angle - start).rem_euclid(TAU);
                if unwrapped > start && unwrapped < end {
                    cuts.push(unwrapped);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        let mut occupied = None;
        for pair in cuts.windows(2) {
            if (pair[1] - pair[0]) * self.radius <= IDENTITY_MM {
                continue;
            }
            let sample = pair[0].midpoint(pair[1]);
            let inside = self
                .arcs
                .iter()
                .filter(|arc| arc.station < level && angle_on_arc(sample, arc.start, arc.end, 0.0))
                .count()
                % 2
                == 1;
            if occupied.is_some_and(|old| old != inside) {
                return None;
            }
            occupied = Some(inside);
        }
        occupied
    }

    const fn curve(&self, start: f64, end: f64) -> RecoveredProfileCurve {
        let (start, end) = if self.outward {
            (start, end)
        } else {
            (end, start)
        };
        RecoveredProfileCurve::CircleArc {
            source_edge_ids: vec![],
            center_mm: self.center,
            radius_mm: self.radius,
            start_angle_rad: start,
            end_angle_rad: end,
        }
    }

    pub(super) fn add_slice(&self, graph: &mut SliceGraph, level: f64) -> Option<()> {
        let mut hits = self
            .generators
            .iter()
            .filter(|g| g.low < level && level < g.high)
            .collect::<Vec<_>>();
        hits.sort_by(|a, b| a.angle.total_cmp(&b.angle));
        if hits.is_empty() {
            if self.occupied(0.0, TAU, level)? {
                graph.closed.push(vec![self.curve(0.0, TAU)]);
            }
            return Some(());
        }
        let mut intervals = Vec::with_capacity(hits.len());
        for (i, hit) in hits.iter().enumerate() {
            let next = hits[(i + 1) % hits.len()];
            let end = if i + 1 == hits.len() {
                next.angle + TAU
            } else {
                next.angle
            };
            if (end - hit.angle) * self.radius <= IDENTITY_MM {
                return None;
            }
            intervals.push((end, self.occupied(hit.angle, end, level)?));
        }
        for (i, hit) in hits.iter().enumerate() {
            let before = intervals[(i + hits.len() - 1) % hits.len()].1;
            let after = intervals[i].1;
            if usize::from(before) + usize::from(after) != hit.multiplicity {
                return None;
            }
            graph.insert(hit.id, point(self.center, self.radius, hit.angle))?;
        }
        for (i, hit) in hits.iter().enumerate() {
            let (end, occupied) = intervals[i];
            if !occupied {
                continue;
            }
            let next = hits[(i + 1) % hits.len()];
            let mut curve = self.curve(hit.angle, end);
            if let RecoveredProfileCurve::CircleArc {
                source_edge_ids, ..
            } = &mut curve
            {
                source_edge_ids.extend([hit.id, next.id]);
                source_edge_ids.sort_unstable();
                source_edge_ids.dedup();
            }
            let (start_id, end_id) = if self.outward {
                (hit.id, next.id)
            } else {
                (next.id, hit.id)
            };
            graph.connect(start_id, end_id, curve)?;
        }
        Some(())
    }
}

fn point(center: [f64; 2], radius: f64, angle: f64) -> [f64; 2] {
    [
        radius.mul_add(angle.cos(), center[0]),
        radius.mul_add(angle.sin(), center[1]),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    fn trim(
        arcs: &[(f64, f64, f64)],
        generators: &[(u64, f64, usize)],
        outward: bool,
    ) -> CylinderTrim {
        CylinderTrim {
            center: [0.0; 2],
            radius: 1.0,
            outward,
            arcs: arcs
                .iter()
                .map(|&(start, end, station)| TransverseArc {
                    center: [0.0; 2],
                    radius: 1.0,
                    start,
                    end,
                    station,
                })
                .collect(),
            generators: generators
                .iter()
                .map(|&(id, angle, multiplicity)| Generator {
                    id,
                    angle,
                    low: 0.0,
                    high: 2.0,
                    multiplicity,
                })
                .collect(),
        }
    }

    #[test]
    fn sectional_cylinder_uses_trim_parity_not_shortest_arc() {
        for (angle, extent) in [(FRAC_PI_2, TAU - FRAC_PI_2), (PI, PI)] {
            let trim = trim(
                &[(angle, TAU, 0.0), (TAU, angle, 2.0)],
                &[(1, 0.0, 1), (2, angle, 1)],
                true,
            );
            let mut graph = SliceGraph::default();
            trim.add_slice(&mut graph, 1.0).unwrap();
            let (end, arc) = &graph.next[&2];
            assert_eq!(*end, 1);
            let RecoveredProfileCurve::CircleArc {
                start_angle_rad,
                end_angle_rad,
                ..
            } = arc
            else {
                panic!("expected circle")
            };
            assert!((end_angle_rad - start_angle_rad - extent).abs() < DIRECTION);
            assert_eq!(graph.next.len(), 1);
        }
    }

    #[test]
    fn sectional_cylinder_periodic_faces_and_seams() {
        for generators in [vec![], vec![(1, 0.0, 2)]] {
            let trim = trim(&[(0.0, TAU, 0.0), (TAU, 0.0, 2.0)], &generators, false);
            let mut graph = SliceGraph::default();
            trim.add_slice(&mut graph, 1.0).unwrap();
            let rings = graph.rings().unwrap();
            assert_eq!(rings.len(), 1);
            assert_eq!(rings[0].len(), 1);
            assert!(crate::solid_extrusions::signed_area_curves(&rings[0]) < 0.0);
        }
    }

    #[test]
    fn sectional_cylinder_rejects_unaccounted_angular_boundary() {
        let trim = trim(&[(0.0, FRAC_PI_2, 0.0), (FRAC_PI_2, 0.0, 2.0)], &[], true);
        let mut graph = SliceGraph::default();
        assert!(trim.add_slice(&mut graph, 1.0).is_none());
    }
}
