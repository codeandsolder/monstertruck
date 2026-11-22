//! Exact correspondence-preserving lofts of planar sections. Repeated
//! polygon points are collapsed slots, not zero-length edges. Constant
//! aligned circular arcs use cylindrical walls; every wall and cap shares
//! the same constructed topology, including triangular limiting walls.
use super::{
    Edge, LoftSection, Point3, Shell, Solid, Wire, builder, curve_end_point, curve_start_point,
    point2_distance, profile, rigid_matrix,
};
use crate::cad_ir::{Curve2d, ProfileLoop};
use anyhow::{Context, Result, ensure};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use monstertruck_modeling::{
    Curve, Face, InnerSpace, Invertible, Line, Matrix4, ParametricSurface3D, Processor,
    RevolutionSurface, Surface, Transform, Vector3, Vertex,
};
use std::{collections::HashMap, f64::consts::TAU};

const JOIN: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.cad_profile_join_mm;

struct Ring {
    vertices: Vec<usize>,
    edges: Vec<Option<Edge>>,
}

type Section = Vec<Ring>;

#[derive(Default)]
struct Skin {
    vertices: Vec<Vertex>,
    edges: HashMap<(usize, usize), Edge>,
}

fn key(p: [f64; 2]) -> [u64; 2] { p.map(|x| if x == 0. { 0 } else { x.to_bits() }) }

/// A collapsed line identifies its two correspondence slots. Prefer the
/// analytic circular endpoint when a limiting corner touches an invariant
/// arc, so both sections use the same transverse cylinder generator.
fn ring_points(curves: &[Curve2d]) -> Result<Vec<[f64; 2]>> {
    let count = curves.len();
    ensure!(count > 0, "empty analytic loft ring");
    let is_arc = |curve: &Curve2d| matches!(curve, Curve2d::CircleArc { .. });
    let candidates = curves
        .iter()
        .enumerate()
        .map(|(i, curve)| {
            let previous = &curves[(i + count - 1) % count];
            let start = curve_start_point(curve)?;
            let end = curve_end_point(previous)?;
            ensure!(
                point2_distance(start, end) <= JOIN,
                "discontinuous analytic loft ring"
            );
            Ok(if !is_arc(curve) && is_arc(previous) {
                (end, true)
            } else {
                (start, is_arc(curve))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let collapsed = |i: usize| match &curves[i] {
        Curve2d::Line { start_mm, end_mm } => key(*start_mm) == key(*end_mm),
        _ => false,
    };
    // Start just after a noncollapsed edge, so a cyclic collapse is
    // represented by one contiguous group rather than two partial groups.
    let offset = (0..count)
        .find(|&i| !collapsed(i))
        .map_or(0, |i| (i + 1) % count);
    let index = |i: usize| (offset + i) % count;
    let mut points = candidates
        .iter()
        .map(|&(point, _)| point)
        .collect::<Vec<_>>();
    let mut first = 0;
    while first < count {
        let mut end = first + 1;
        while end < count && collapsed(index(end - 1)) {
            end += 1;
        }
        let anchor = (first..end)
            .map(index)
            .find(|&i| candidates[i].1)
            .unwrap_or_else(|| index(first));
        let canonical = candidates[anchor].0;
        for i in (first..end).map(index) {
            ensure!(
                point2_distance(points[i], canonical) <= JOIN,
                "inconsistent collapsed loft corner"
            );
            points[i] = canonical;
        }
        first = end;
    }
    for (i, curve) in curves.iter().enumerate() {
        ensure!(
            point2_distance(curve_start_point(curve)?, points[i]) <= JOIN
                && point2_distance(curve_end_point(curve)?, points[(i + 1) % count]) <= JOIN,
            "canonical loft vertex outside its source curve endpoints"
        );
    }
    Ok(points)
}

impl Skin {
    fn section(&mut self, section: &LoftSection, profiles: &[ProfileLoop]) -> Result<Section> {
        let matrix = rigid_matrix(&section.transform);
        profiles
            .iter()
            .map(|ring| self.ring(ring, matrix))
            .collect()
    }

    fn ring(&mut self, ring: &ProfileLoop, matrix: Matrix4) -> Result<Ring> {
        let curves = &ring.curves;
        let points = ring_points(curves)?;
        let mut indices = HashMap::new();
        let vertices = points
            .iter()
            .map(|&p| {
                *indices.entry(key(p)).or_insert_with(|| {
                    let index = self.vertices.len();
                    self.vertices.push(builder::vertex(
                        matrix.transform_point(Point3::new(p[0], p[1], 0.)),
                    ));
                    index
                })
            })
            .collect::<Vec<_>>();
        let edges = curves
            .iter()
            .enumerate()
            .map(|(i, curve)| {
                let a = vertices[i];
                let b = vertices[(i + 1) % vertices.len()];
                match curve {
                    Curve2d::Line { .. } => Ok((a != b).then(|| self.edge(a, b))),
                    Curve2d::CircleArc {
                        center_mm,
                        radius_mm,
                        start_angle_rad,
                        end_angle_rad,
                    } => {
                        ensure!(a != b, "collapsed circular loft edge");
                        let mid = start_angle_rad.midpoint(*end_angle_rad);
                        let transit = matrix.transform_point(Point3::new(
                            radius_mm.mul_add(mid.cos(), center_mm[0]),
                            radius_mm.mul_add(mid.sin(), center_mm[1]),
                            0.,
                        ));
                        Ok(Some(builder::circle_arc(
                            &self.vertices[a],
                            &self.vertices[b],
                            transit,
                        )?))
                    }
                    _ => anyhow::bail!("unsupported analytic loft curve"),
                }
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Ring { vertices, edges })
    }

    fn edge(&mut self, a: usize, b: usize) -> Edge {
        let (key, forward) = if a < b {
            ((a, b), true)
        } else {
            ((b, a), false)
        };
        let edge = self
            .edges
            .entry(key)
            .or_insert_with(|| builder::line(&self.vertices[key.0], &self.vertices[key.1]));
        if forward {
            edge.clone()
        } else {
            edge.inverse()
        }
    }

    fn cap(shell: &mut Shell, section: &Section, reverse: bool) -> Result<()> {
        let wires = section
            .iter()
            .filter_map(|ring| {
                let edges = ring.edges.iter().flatten().cloned().collect::<Vec<_>>();
                // A polygon can collapse to a line or point. A two-edge
                // circular segment, however, has positive area.
                (edges.len() >= 3
                    || (edges.len() == 2
                        && edges
                            .iter()
                            .any(|edge| !matches!(edge.curve(), Curve::Line(_)))))
                .then(|| Wire::from(edges))
            })
            .collect::<Vec<_>>();
        if wires.is_empty() {
            return Ok(());
        }
        let cap = profile::attach_plane_normalized(wires)?;
        shell.push(if reverse { cap.inverse() } else { cap });
        Ok(())
    }

    fn connect(
        &mut self,
        shell: &mut Shell,
        start: &Section,
        end: &Section,
        section: &LoftSection,
        profiles: &[ProfileLoop],
    ) -> Result<()> {
        ensure!(
            start.len() == end.len() && start.len() == profiles.len(),
            "incompatible analytic loft ring count"
        );
        for ((a, b), profile) in start.iter().zip(end).zip(profiles) {
            ensure!(
                [a.vertices.len(), b.vertices.len()]
                    .into_iter()
                    .all(|len| len == profile.curves.len()),
                "incompatible analytic loft correspondence"
            );
            for (i, curve) in profile.curves.iter().enumerate() {
                let j = (i + 1) % a.vertices.len();
                if a.edges[i].is_none() && b.edges[i].is_none() {
                    continue;
                }
                let mut edges = Vec::with_capacity(4);
                edges.extend(a.edges[i].clone());
                edges.push(self.edge(a.vertices[j], b.vertices[j]));
                edges.extend(b.edges[i].as_ref().map(Edge::inverse));
                edges.push(self.edge(b.vertices[i], a.vertices[i]));
                let wire = Wire::from(edges);
                let face =
                    match curve {
                        Curve2d::Line { .. } => builder::try_attach_plane(vec![wire])?,
                        Curve2d::CircleArc {
                            center_mm,
                            start_angle_rad,
                            end_angle_rad,
                            ..
                        } => {
                            let low = self.vertices[a.vertices[i]].point();
                            let high = self.vertices[b.vertices[i]].point();
                            let axis = (high - low).normalize();
                            let center = rigid_matrix(&section.transform)
                                .transform_point(Point3::new(center_mm[0], center_mm[1], 0.));
                            let radial = low - center;
                            let forward = normal(section).dot(axis) > 0.;
                            let outward = (*end_angle_rad > *start_angle_rad) == forward;
                            let mut surface = Surface::Revolution(Processor::new(
                                RevolutionSurface::by_revolution(
                                    Curve::Line(Line(low, high)),
                                    center,
                                    axis,
                                ),
                            ));
                            if (surface.normal(0.5, 0.).dot(radial) > 0.) != outward {
                                surface.invert();
                            }
                            Face::new(vec![wire], surface)?
                        }
                        _ => anyhow::bail!("unsupported analytic loft wall"),
                    };
                shell.push(face);
            }
        }
        Ok(())
    }
}

const fn normal(s: &LoftSection) -> Vector3 {
    Vector3::new(
        s.transform.matrix[0][2],
        s.transform.matrix[1][2],
        s.transform.matrix[2][2],
    )
}
const fn origin(s: &LoftSection) -> Vector3 {
    Vector3::new(
        s.transform.matrix[0][3],
        s.transform.matrix[1][3],
        s.transform.matrix[2][3],
    )
}

/// Signed area in the section chart, including the exact circular segment
/// area rather than a chord-only winding approximation.
fn signed_area(ring: &ProfileLoop) -> Result<f64> {
    let first = curve_start_point(ring.curves.first().context("empty loft ring")?)?;
    ring.curves
        .iter()
        .map(|curve| {
            let a = curve_start_point(curve)?;
            let b = curve_end_point(curve)?;
            let [ax, ay] = [a[0] - first[0], a[1] - first[1]];
            let [bx, by] = [b[0] - first[0], b[1] - first[1]];
            let segment = match curve {
                Curve2d::CircleArc {
                    radius_mm,
                    start_angle_rad,
                    end_angle_rad,
                    ..
                } => {
                    let sweep = end_angle_rad - start_angle_rad;
                    radius_mm.powi(2) * (sweep - sweep.sin())
                }
                _ => 0.,
            };
            Ok(ax.mul_add(by, -ay * bx) + segment)
        })
        .sum::<Result<f64>>()
}

fn reverse_ring(ring: &mut ProfileLoop) {
    ring.curves.reverse();
    for curve in &mut ring.curves {
        match curve {
            Curve2d::Line { start_mm, end_mm } => core::mem::swap(start_mm, end_mm),
            Curve2d::CircleArc {
                start_angle_rad,
                end_angle_rad,
                ..
            } => core::mem::swap(start_angle_rad, end_angle_rad),
            _ => {}
        }
    }
}

fn orient_profiles(sections: &[LoftSection], profiles: &mut [Vec<ProfileLoop>]) -> Result<()> {
    let axis = normal(&sections[0]);
    let progression =
        origin(sections.last().context("missing final section")?) - origin(&sections[0]);
    let forward = axis.dot(progression) > 0.;
    let areas = profiles
        .iter()
        .map(|rings| rings.iter().map(signed_area).collect::<Result<Vec<_>>>())
        .collect::<Result<Vec<_>>>()?;
    let first = areas
        .iter()
        .find_map(|p| (p[0] != 0.).then_some(p[0] > 0.))
        .context("all loft outer sections have zero area")?;
    let reverse = first != forward;
    for (i, section) in sections.iter().enumerate() {
        ensure!(
            normal(section).dot(axis) > 0.
                && (i == 0
                    || (origin(section) - origin(&sections[i - 1])).dot(axis)
                        * (if forward { 1. } else { -1. })
                        > 0.),
            "analytic loft sections lack consistently oriented axial progression"
        );
        for (j, (ring, &area)) in profiles[i].iter_mut().zip(&areas[i]).enumerate() {
            ensure!(
                area == 0. || (area > 0.) == (if j == 0 { first } else { !first }),
                "inconsistent loft ring winding"
            );
            if reverse {
                reverse_ring(ring);
            }
        }
    }
    Ok(())
}

/// Full circles use two correspondence slots, keeping the cap and
/// cylindrical skin topology ordinary and nondegenerate.
fn split_circles(ring: &ProfileLoop) -> ProfileLoop {
    ProfileLoop {
        curves: ring
            .curves
            .iter()
            .flat_map(|curve| {
                if let Curve2d::CircleArc {
                    center_mm,
                    radius_mm,
                    start_angle_rad,
                    end_angle_rad,
                } = curve
                    && ((end_angle_rad - start_angle_rad).abs() - TAU).abs()
                        <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
                {
                    let mid = start_angle_rad.midpoint(*end_angle_rad);
                    vec![
                        Curve2d::CircleArc {
                            center_mm: *center_mm,
                            radius_mm: *radius_mm,
                            start_angle_rad: *start_angle_rad,
                            end_angle_rad: mid,
                        },
                        Curve2d::CircleArc {
                            center_mm: *center_mm,
                            radius_mm: *radius_mm,
                            start_angle_rad: mid,
                            end_angle_rad: *end_angle_rad,
                        },
                    ]
                } else {
                    vec![curve.clone()]
                }
            })
            .collect(),
    }
}

fn aligned_circular_sections(sections: &[LoftSection], profiles: &[Vec<ProfileLoop>]) -> bool {
    sections.windows(2).zip(profiles.windows(2)).all(|(s, p)| {
        let axis = normal(&s[0]);
        let delta = origin(&s[1]) - origin(&s[0]);
        // Only unchanged local arcs with a pure axial translation define
        // an exact cylinder. A twisted or varying-radius loft belongs to
        // the general homotopy skin instead.
        (0..3).all(|i| {
            (0..3).all(|j| {
                s[0].transform.matrix[i][j].to_bits() == s[1].transform.matrix[i][j].to_bits()
            })
        }) && delta.cross(axis).magnitude() <= JOIN
            && p[0].len() == p[1].len()
            && p[0].iter().zip(&p[1]).all(|(a, b)| {
                a.curves.len() == b.curves.len()
                    && a.curves.iter().zip(&b.curves).all(|(a, b)| match (a, b) {
                        (Curve2d::Line { .. }, Curve2d::Line { .. }) => true,
                        (Curve2d::CircleArc { .. }, Curve2d::CircleArc { .. }) => a == b,
                        _ => false,
                    })
            })
    })
}

pub(super) fn build_if_needed(sections: &[LoftSection]) -> Result<Option<Solid>> {
    let profiles = sections.iter().map(|s| &s.profile.loops);
    let all_lines = profiles
        .clone()
        .flatten()
        .all(|ring| ring.polygon_points().is_some());
    let supported = profiles
        .clone()
        .flatten()
        .flat_map(|ring| &ring.curves)
        .all(|curve| matches!(curve, Curve2d::Line { .. } | Curve2d::CircleArc { .. }));
    if !supported || profiles.clone().any(Vec::is_empty) {
        return Ok(None);
    }
    let needed = profiles.clone().any(|rings| {
        rings.len() > 1
            || rings
                .iter()
                .flat_map(|ring| &ring.curves)
                .any(|curve| match curve {
                    Curve2d::Line { start_mm, end_mm } => key(*start_mm) == key(*end_mm),
                    _ => false,
                })
    });
    if !needed {
        return Ok(None);
    }
    let mut profiles = profiles
        .map(|rings| rings.iter().map(split_circles).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    if !all_lines && !aligned_circular_sections(sections, &profiles) {
        return Ok(None);
    }
    ensure!(
        profiles
            .windows(2)
            .all(|pair| pair[0].len() == pair[1].len()
                && pair[0]
                    .iter()
                    .zip(&pair[1])
                    .all(|(a, b)| a.curves.len() == b.curves.len())),
        "incompatible analytic loft profiles"
    );
    orient_profiles(sections, &mut profiles)?;
    let mut skin = Skin::default();
    let topology = sections
        .iter()
        .zip(&profiles)
        .map(|(s, p)| skin.section(s, p))
        .collect::<Result<Vec<_>>>()?;
    let mut shell = Shell::new();
    Skin::cap(&mut shell, &topology[0], true)?;
    for (i, pair) in topology.windows(2).enumerate() {
        skin.connect(&mut shell, &pair[0], &pair[1], &sections[i], &profiles[i])?;
    }
    Skin::cap(
        &mut shell,
        topology.last().context("missing final topology")?,
        false,
    )?;
    Ok(Some(Solid::new(vec![shell])?))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cad_ir::{Profile2d, RigidTransform};
    use monstertruck_modeling::MetricSpace;

    const fn section(profile: Profile2d, z: f64) -> LoftSection {
        LoftSection {
            profile,
            transform: RigidTransform::translation_mm([0.0, 0.0, z]),
        }
    }

    #[test]
    fn limiting_planar_loft_builds_triangular_walls() -> Result<()> {
        let base = Profile2d::polygon(vec![[0., 0.], [4., 0.], [3., 2.], [1., 2.]])?;
        let tip = Profile2d::polygon(vec![[0., 0.], [4., 0.], [2., 4.], [2., 4.]])?;
        let solid = build_if_needed(&[section(base, 0.0), section(tip, 2.0)])?
            .context("limiting loft expected")?;
        let shells = solid.boundaries();
        let [shell] = shells.as_slice() else {
            anyhow::bail!("expected one shell");
        };
        assert!(solid.is_geometric_consistent());
        assert_eq!(shell.len(), 6);
        assert_eq!(
            shell
                .iter()
                .filter(|f| f.boundaries()[0].len() == 3)
                .count(),
            2
        );
        assert!(
            shell.iter().flat_map(Face::edge_iter).all(|e| e
                .front()
                .point()
                .distance(e.back().point())
                > 0.0)
        );
        Ok(())
    }

    #[test]
    fn limiting_loft_preserves_circular_hole_and_collapsed_planar_corner() -> Result<()> {
        let hollow = |points: Vec<[f64; 2]>| -> Result<Profile2d> {
            let mut profile = Profile2d::polygon(points)?;
            profile.loops.push(ProfileLoop {
                curves: vec![Curve2d::CircleArc {
                    center_mm: [2., 1.],
                    radius_mm: 0.3,
                    start_angle_rad: 0.,
                    end_angle_rad: -TAU,
                }],
            });
            Ok(profile)
        };
        let base = hollow(vec![[0., 0.], [4., 0.], [3., 2.], [1., 2.]])?;
        let tip = hollow(vec![[0., 0.], [4., 0.], [2., 4.], [2., 4.]])?;
        let solid = build_if_needed(&[section(base, 0.), section(tip, 2.)])?
            .context("mixed limiting loft expected")?;
        let shells = solid.boundaries();
        let [shell] = shells.as_slice() else {
            anyhow::bail!("expected one shell");
        };
        assert!(solid.is_geometric_consistent());
        assert_eq!(
            shell
                .iter()
                .filter(|face| matches!(face.surface(), Surface::Revolution(_)))
                .count(),
            2
        );
        assert_eq!(
            shell
                .iter()
                .filter(|face| face.boundaries().len() == 2)
                .count(),
            2
        );
        assert!(
            shell.iter().flat_map(Face::edge_iter).all(|edge| edge
                .front()
                .point()
                .distance(edge.back().point())
                > 0.)
        );
        Ok(())
    }

    #[test]
    fn planar_loft_keeps_tapered_hole_in_both_caps() -> Result<()> {
        let base = Profile2d::polygon(vec![[-3., -3.], [3., -3.], [3., 3.], [-3., 3.]])?;
        let hollow = |r: f64| -> Result<Profile2d> {
            let mut profile = base.clone();
            profile
                .loops
                .extend(Profile2d::polygon(vec![[-r, -r], [-r, r], [r, r], [r, -r]])?.loops);
            Ok(profile)
        };
        let solid = build_if_needed(&[section(hollow(1.0)?, 0.), section(hollow(0.5)?, 2.)])?
            .context("hollow loft expected")?;
        let shells = solid.boundaries();
        let [shell] = shells.as_slice() else {
            anyhow::bail!("expected single connected shell");
        };
        assert!(solid.is_geometric_consistent());
        assert_eq!(shell.len(), 10);
        assert_eq!(
            shell.iter().filter(|f| f.boundaries().len() == 2).count(),
            2
        );
        Ok(())
    }
}
