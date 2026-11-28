//! Analytic swept slabs with linearly changing planar edges and fixed arcs. Within consecutive
//! vertex stations, each finite face intersects a slice in the same ordered
//! intervals. The mid-slice graph supplies correspondence; the actual
//! original crossing edges supply both end profiles, including collapsed
//! corner transitions. Axis-parallel cylinders have invariant section arcs. Independent full-source proof remains mandatory.
use super::{
    DIRECTION, Frame, IDENTITY_MM, circular::CylinderTrim, section_groups, slice_rings,
    station_levels, valid_planar_face,
};
use crate::{
    brep::{OrientedEdgeUse, SurfaceSupport},
    cad_ir::{CadModel, Curve2d, NodeId, Profile2d, ProfileLoop},
    cad_recovery::{add_bounded_union_tree, solid::recovered_profile_curve_to_ir},
    math2::distance as distance2,
    math3::{add, cross, dot, mul, norm, normalize, sub},
    profile_curves::RecoveredProfileCurve,
    reduced_brep::{ReducedFace, ReducedSolid},
};
use anyhow::Result;
use std::collections::BTreeMap;

mod chains;

struct RuledSweep<'a> {
    edges: BTreeMap<u64, &'a OrientedEdgeUse>,
    levels: Vec<f64>,
}

impl<'a> RuledSweep<'a> {
    fn new(solid: &'a ReducedSolid, frame: Frame) -> Option<Self> {
        let mut inclined = false;
        for face in solid.faces() {
            match face.surface {
                SurfaceSupport::Plane(plane) => {
                    let normal = normalize(plane.normal, DIRECTION)?;
                    inclined |= norm(cross(normal, frame.axis)) > DIRECTION
                        && dot(normal, frame.axis).abs() > DIRECTION;
                    valid_planar_face(face, frame, true)?;
                }
                SurfaceSupport::Cylinder(_) => {
                    CylinderTrim::new(face, frame)?;
                }
                _ => return None,
            }
        }
        if !inclined {
            return None;
        }
        let levels = station_levels(solid, frame)?;
        let edges = solid
            .faces()
            .flat_map(ReducedFace::loop_edges)
            .map(|edge| (edge.edge_id, edge))
            .collect();
        Some(Self { edges, levels })
    }

    fn point(&self, id: u64, frame: Frame, level: f64) -> Option<[f64; 2]> {
        let edge = self.edges.get(&id)?;
        let [a, b] = [edge.start_mm, edge.end_mm].map(|p| frame.axial(p));
        if (a - b).abs() <= IDENTITY_MM
            || level < a.min(b) - IDENTITY_MM
            || level > a.max(b) + IDENTITY_MM
        {
            return None;
        }
        let point = if (a - level).abs() <= IDENTITY_MM {
            edge.start_mm
        } else if (b - level).abs() <= IDENTITY_MM {
            edge.end_mm
        } else {
            add(
                edge.start_mm,
                mul(sub(edge.end_mm, edge.start_mm), (level - a) / (b - a)),
            )
        };
        Some(frame.project(point))
    }

    fn profile(&self, rings: &Section, frame: Frame, level: f64) -> Option<Profile2d> {
        let loops = rings
            .iter()
            .map(|ring| {
                let mut points = ring
                    .iter()
                    .map(|curve| match curve {
                        RecoveredProfileCurve::Line {
                            source_edge_ids, ..
                        } => {
                            let [start, _] = source_edge_ids.as_slice() else {
                                return None;
                            };
                            self.point(*start, frame, level)
                        }
                        RecoveredProfileCurve::CircleArc { .. } => curve.start_point(),
                        _ => None,
                    })
                    .collect::<Option<Vec<_>>>()?;
                if points.len() >= 3 {
                    canonicalize_collapsed_corners(&mut points)?;
                }
                let curves = ring
                    .iter()
                    .zip(points.iter().zip(points.iter().cycle().skip(1)))
                    .map(|(curve, (&start_mm, &end_mm))| match curve {
                        RecoveredProfileCurve::Line {
                            source_edge_ids, ..
                        } => {
                            let [_, end] = source_edge_ids.as_slice() else {
                                return None;
                            };
                            (distance2(end_mm, self.point(*end, frame, level)?) <= IDENTITY_MM)
                                .then_some(Curve2d::Line { start_mm, end_mm })
                        }
                        RecoveredProfileCurve::CircleArc { .. } => {
                            (distance2(start_mm, curve.start_point()?) <= IDENTITY_MM
                                && distance2(end_mm, curve.end_point()?) <= IDENTITY_MM)
                                .then(|| recovered_profile_curve_to_ir(curve))
                        }
                        _ => None,
                    })
                    .collect::<Option<Vec<_>>>()?;
                Some(ProfileLoop { curves })
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Profile2d { loops })
    }
}

/// Canonicalize contiguous limiting corners, using the original point
/// distances rather than the progressively modified coordinates. Every
/// displacement must fit the existing length budget, even along a chain
/// of individually close points.
fn canonicalize_collapsed_corners(points: &mut [[f64; 2]]) -> Option<()> {
    let count = points.len();
    if count < 3 {
        return None;
    }
    let offset = (0..count)
        .find(|&i| distance2(points[i], points[(i + 1) % count]) > IDENTITY_MM)
        .map_or(0, |i| (i + 1) % count);
    points.rotate_left(offset);
    let mut start = 0;
    while start < count {
        let mut end = start + 1;
        while end < count && distance2(points[end - 1], points[end]) <= IDENTITY_MM {
            end += 1;
        }
        let group = &mut points[start..end];
        let canonical = *group
            .iter()
            .min_by(|a, b| a[0].total_cmp(&b[0]).then_with(|| a[1].total_cmp(&b[1])))?;
        if group.iter().any(|&p| distance2(p, canonical) > IDENTITY_MM) {
            return None;
        }
        group.fill(canonical);
        start = end;
    }
    points.rotate_right(offset);
    Some(())
}

/// Complete ordered crossing-edge correspondence of an outer ring and its
/// holes. An unchanged key across source stations establishes identical
/// affine generators throughout the interval; no profile fitting or
/// near-coordinate matching is needed.
type Section = Vec<Vec<RecoveredProfileCurve>>;
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum CurveKey {
    Line([u64; 2]),
    Arc {
        source_edges: Vec<u64>,
        center: [u64; 2],
        radius: u64,
        angles: [u64; 2],
    },
}

impl CurveKey {
    fn new(curve: &RecoveredProfileCurve) -> Option<Self> {
        match curve {
            RecoveredProfileCurve::Line {
                source_edge_ids, ..
            } => {
                let [start, end] = source_edge_ids.as_slice() else {
                    return None;
                };
                Some(Self::Line([*start, *end]))
            }
            RecoveredProfileCurve::CircleArc {
                source_edge_ids,
                center_mm,
                radius_mm,
                start_angle_rad,
                end_angle_rad,
            } => Some(Self::Arc {
                source_edges: source_edge_ids.clone(),
                center: center_mm.map(f64::to_bits),
                radius: radius_mm.to_bits(),
                angles: [start_angle_rad.to_bits(), end_angle_rad.to_bits()],
            }),
            _ => None,
        }
    }
}
type SectionKey = Vec<Vec<CurveKey>>;

fn slab_sections(
    solid: &ReducedSolid,
    frame: Frame,
    limits: [f64; 2],
) -> Option<Vec<(SectionKey, Section)>> {
    let rings = slice_rings(solid, frame, limits[0].midpoint(limits[1]))?;
    if rings.is_empty() {
        return Some(Vec::new());
    }
    section_groups(&rings)?
        .iter()
        .map(|group| {
            let section: Section = group.iter().map(|&i| rings[i].clone()).collect();
            let key = section
                .iter()
                .map(|ring| ring.iter().map(CurveKey::new).collect())
                .collect::<Option<_>>()?;
            Some((key, section))
        })
        .collect()
}

/// Maintain maximal contiguous station intervals for each unchanged
/// connected section. A different island can start, stop or change shape
/// without splitting unaffected material into additional Boolean operands.
#[derive(Default)]
struct SectionRuns {
    active: BTreeMap<SectionKey, usize>,
    completed: Vec<(usize, usize, SectionKey)>,
}

impl SectionRuns {
    fn advance(&mut self, station: usize, sections: Vec<SectionKey>) {
        let next = sections
            .into_iter()
            .map(|section| {
                let start = self.active.remove(&section).unwrap_or(station);
                (section, start)
            })
            .collect();
        self.completed.extend(
            std::mem::take(&mut self.active)
                .into_iter()
                .map(|(section, start)| (start, station, section)),
        );
        self.active = next;
    }

    fn finish(mut self, station: usize) -> Vec<(usize, usize, SectionKey)> {
        self.advance(station, Vec::new());
        self.completed.sort_unstable();
        self.completed
    }
}

fn same_profiles(a: &Profile2d, b: &Profile2d) -> bool {
    a.loops.len() == b.loops.len() && a.loops.iter().zip(&b.loops).all(|(a, b)| {
        a.curves.len() == b.curves.len() && a.curves.iter().zip(&b.curves).all(|(a,b)| {
            !matches!(a, Curve2d::Line { start_mm, end_mm } if distance2(*start_mm,*end_mm) <= IDENTITY_MM)
                && chains::matching_curve(a,b)
        })
    })
}

fn section_spans(
    sweep: &RuledSweep<'_>,
    frame: Frame,
    intervals: Vec<(usize, usize, SectionKey)>,
    sections: &BTreeMap<SectionKey, Section>,
) -> Option<Vec<chains::Span>> {
    intervals
        .into_iter()
        .map(|(lo, hi, section)| {
            Some(chains::Span {
                limits: [lo, hi],
                start: sweep.profile(sections.get(&section)?, frame, sweep.levels[lo])?,
                end: sweep.profile(sections.get(&section)?, frame, sweep.levels[hi])?,
            })
        })
        .collect()
}

fn emit_model(
    chains: Vec<chains::Chain>,
    frame: Frame,
    levels: &[f64],
) -> Result<Option<(CadModel, NodeId, usize)>> {
    let mut model = CadModel::new();
    let mut grouped = BTreeMap::<usize, Vec<NodeId>>::new();
    for chain in chains {
        let lo = chain.first_station();
        let node = chain.into_node(&mut model, frame, levels);
        grouped.entry(lo).or_default().push(node);
    }
    let children = grouped
        .into_values()
        .map(|parts| add_bounded_union_tree(&mut model, parts))
        .collect::<Vec<_>>();
    if children.is_empty() {
        return Ok(None);
    }
    let root = add_bounded_union_tree(&mut model, children);
    model.add_root(root)?;
    model.validate()?;
    Ok(Some((model, root, levels.len() - 1)))
}

/// Keep the original slab construction eligible for selection. The
/// independently certified candidate ranking, not a heuristic preference
/// for loft chains, decides whether a compact alternative wins.
pub(super) fn models(solid: &ReducedSolid, frame: Frame) -> Result<Vec<(CadModel, NodeId, usize)>> {
    let Some(sweep) = RuledSweep::new(solid, frame) else {
        return Ok(Vec::new());
    };
    let mut runs = SectionRuns::default();
    let mut slabs = Vec::new();
    let mut section_table = BTreeMap::new();
    for (station, ab) in sweep.levels.windows(2).enumerate() {
        let Some(regions) = slab_sections(solid, frame, [ab[0], ab[1]]) else {
            return Ok(Vec::new());
        };
        let sections = regions
            .into_iter()
            .map(|(key, section)| {
                section_table.entry(key.clone()).or_insert(section);
                key
            })
            .collect::<Vec<_>>();
        slabs.extend(
            sections
                .iter()
                .cloned()
                .map(|section| (station, station + 1, section)),
        );
        runs.advance(station, sections);
    }
    let coalesced = runs.finish(sweep.levels.len() - 1);
    let distinct = coalesced.len() < slabs.len();
    let Some(baseline) = section_spans(&sweep, frame, slabs, &section_table) else {
        return Ok(Vec::new());
    };
    let mut models = Vec::new();
    if let Some(model) = emit_model(
        baseline.iter().cloned().map(chains::Chain::new).collect(),
        frame,
        &sweep.levels,
    )? {
        models.push(model);
    }
    let spans = if distinct {
        let Some(spans) = section_spans(&sweep, frame, coalesced, &section_table) else {
            return Ok(models);
        };
        if let Some(model) = emit_model(
            spans.iter().cloned().map(chains::Chain::new).collect(),
            frame,
            &sweep.levels,
        )? {
            models.push(model);
        }
        spans
    } else {
        baseline
    };
    let count = spans.len();
    let joined = chains::join(spans);
    if joined.len() < count
        && let Some(model) = emit_model(joined, frame, &sweep.levels)?
    {
        models.push(model);
    }
    Ok(models)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unrelated_stations_do_not_split_unchanged_sections() {
        let first = vec![vec![
            CurveKey::Line([1, 2]),
            CurveKey::Line([2, 3]),
            CurveKey::Line([3, 4]),
            CurveKey::Line([4, 1]),
        ]];
        let second = vec![vec![
            CurveKey::Line([5, 6]),
            CurveKey::Line([6, 7]),
            CurveKey::Line([7, 5]),
        ]];
        let mut runs = SectionRuns::default();
        runs.advance(0, vec![first.clone()]);
        runs.advance(1, vec![first.clone(), second.clone()]);
        runs.advance(2, vec![second.clone()]);
        runs.advance(3, vec![first.clone(), second.clone()]);
        assert_eq!(
            runs.finish(4),
            vec![(0, 2, first.clone()), (1, 4, second), (3, 4, first)]
        );
    }

    #[test]
    fn collapsed_corners_preserve_cyclic_correspondence() -> Result<()> {
        let mut points = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [IDENTITY_MM * 0.2, 0.0]];
        anyhow::Context::context(
            canonicalize_collapsed_corners(&mut points),
            "valid cyclic collapse",
        )?;
        assert!(distance2(points[0], points[3]) <= f64::EPSILON);
        assert!(distance2(points[1], [2.0, 0.0]) <= f64::EPSILON);
        assert!(distance2(points[2], [2.0, 2.0]) <= f64::EPSILON);
        Ok(())
    }

    #[test]
    fn collapsed_corners_reject_transitive_displacement() {
        let mut points = [
            [0.0, 0.0],
            [0.75 * IDENTITY_MM, 0.0],
            [1.5 * IDENTITY_MM, 0.0],
            [2.0, 2.0],
        ];
        assert!(canonicalize_collapsed_corners(&mut points).is_none());
    }
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod native_tests {
    use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
    use crate::complete_ir::recover_complete_ir_with_kernel_bytes;
    use anyhow::Result;

    #[test]
    fn tapered_host_prefers_a_certified_multistation_loft() -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let bytes =
            include_bytes!("../../../validation/fixtures/native_peeled_tapered_host_6629a4db.step");
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        let fragments = report.selected_fragments().collect::<Vec<_>>();
        assert_eq!(fragments.len(), 1);
        let fragment = &fragments[0].fragment;
        // Measured original independently certified slab cost: 526.
        assert!(fragment.model.complexity_score(fragment.root)? < 526);
        assert!(fragment.model.nodes.iter().any(
            |node| matches!(node, crate::cad_ir::CadNode::Loft { sections } if sections.len() > 2)
        ));
        Ok(())
    }

    #[test]
    fn inclined_relief_host_recovery_is_source_verified() -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let bytes =
            include_bytes!("../../../validation/fixtures/native_relief_oblique_8ddfec04.step");
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        let fragments = report.selected_fragments().collect::<Vec<_>>();
        assert_eq!(fragments.len(), 1);
        let fragment = &fragments[0].fragment;
        let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
        let summary = kernel.summarize(&evaluated);
        assert!(summary.geometrically_consistent);
        assert_eq!(summary.shells, 1);
        Ok(())
    }
}
