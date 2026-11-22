//! Exact correspondence-preserving planar loft skins. Repeated polygon
//! points represent collapsed correspondence slots, not zero-length edges.
//! Construct triangular limiting walls and reuse every actual topological
//! edge between walls and caps.
use super::{Edge, LoftSection, Point3, Shell, Solid, Wire, builder, profile, rigid_matrix};
use crate::cad_ir::ProfileLoop;
use anyhow::{Context, Result, ensure};
use monstertruck_modeling::{InnerSpace, Vector3, Vertex};
use std::collections::HashMap;

type Ring = Vec<usize>;
type Section = Vec<Ring>;

#[derive(Default)]
struct Skin {
    vertices: Vec<Vertex>,
    edges: HashMap<(usize, usize), Edge>,
}

fn key(p: [f64; 2]) -> [u64; 2] { p.map(|x| if x == 0. { 0 } else { x.to_bits() }) }

impl Skin {
    fn section(&mut self, section: &LoftSection, polygons: &[Vec<[f64; 2]>]) -> Section {
        let matrix = rigid_matrix(&section.transform);
        polygons
            .iter()
            .map(|ring| {
                let mut indices = HashMap::new();
                ring.iter()
                    .map(|&p| {
                        *indices.entry(key(p)).or_insert_with(|| {
                            let index = self.vertices.len();
                            self.vertices.push(builder::transformed(
                                &builder::vertex(Point3::new(p[0], p[1], 0.)),
                                matrix,
                            ));
                            index
                        })
                    })
                    .collect()
            })
            .collect()
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

    fn wire(&mut self, ring: &[usize]) -> Option<Wire> {
        let mut vertices = ring.to_vec();
        vertices.dedup();
        if vertices.first() == vertices.last() {
            vertices.pop();
        }
        if vertices.len() < 3 {
            return None;
        }
        let edges = vertices
            .iter()
            .zip(vertices.iter().cycle().skip(1))
            .map(|(&a, &b)| self.edge(a, b))
            .collect::<Vec<_>>();
        Some(Wire::from(edges))
    }

    fn cap(&mut self, shell: &mut Shell, section: &Section, reverse: bool) -> Result<()> {
        let wires = section
            .iter()
            .filter_map(|ring| self.wire(ring))
            .collect::<Vec<_>>();
        if wires.is_empty() {
            return Ok(());
        }
        let cap = profile::attach_plane_normalized(wires)?;
        shell.push(if reverse { cap.inverse() } else { cap });
        Ok(())
    }

    fn connect(&mut self, shell: &mut Shell, start: &Section, end: &Section) -> Result<()> {
        for (a, b) in start.iter().zip(end) {
            ensure!(
                a.len() == b.len(),
                "incompatible planar loft correspondence"
            );
            for i in 0..a.len() {
                let j = (i + 1) % a.len();
                if let Some(wire) = self.wire(&[a[i], a[j], b[j], b[i]]) {
                    shell.push(builder::try_attach_plane(vec![wire])?);
                }
            }
        }
        Ok(())
    }
}

fn signed_area(p: &[[f64; 2]]) -> f64 {
    // Shift before the cross products to avoid cancellation for far-away
    // placements. This establishes winding, not a geometric fit.
    let origin = p[0];
    p.iter()
        .zip(p.iter().cycle().skip(1))
        .map(|(a, b)| {
            let [ax, ay] = [a[0] - origin[0], a[1] - origin[1]];
            let [bx, by] = [b[0] - origin[0], b[1] - origin[1]];
            ax.mul_add(by, -ay * bx)
        })
        .sum::<f64>()
}

fn orient_polygons(sections: &[LoftSection], polygons: &mut [Vec<Vec<[f64; 2]>>]) -> Result<()> {
    let normal = |s: &LoftSection| {
        Vector3::new(
            s.transform.matrix[0][2],
            s.transform.matrix[1][2],
            s.transform.matrix[2][2],
        )
    };
    let origin = |s: &LoftSection| {
        Vector3::new(
            s.transform.matrix[0][3],
            s.transform.matrix[1][3],
            s.transform.matrix[2][3],
        )
    };
    let axis = normal(&sections[0]);
    let progression =
        origin(sections.last().context("missing final section")?) - origin(&sections[0]);
    let forward = axis.dot(progression) > 0.;
    let first = polygons
        .iter()
        .find_map(|p| {
            let area = signed_area(&p[0]);
            (area != 0.).then_some(area > 0.)
        })
        .context("all loft outer sections have zero area")?;
    let reverse = first != forward;
    for (i, section) in sections.iter().enumerate() {
        ensure!(
            normal(section).dot(axis) > 0.
                && (i == 0
                    || (origin(section) - origin(&sections[i - 1])).dot(axis)
                        * (if forward { 1. } else { -1. })
                        > 0.),
            "planar loft sections do not have consistently oriented axial progression"
        );
        for (j, ring) in polygons[i].iter_mut().enumerate() {
            let area = signed_area(ring);
            if area != 0. {
                ensure!(
                    (area > 0.) == (if j == 0 { first } else { !first }),
                    "inconsistent loft ring winding"
                );
            }
            if reverse {
                ring.reverse();
            }
        }
    }
    Ok(())
}

pub(super) fn build_if_needed(sections: &[LoftSection]) -> Result<Option<Solid>> {
    let Some(mut polygons) = sections
        .iter()
        .map(|s| {
            s.profile
                .loops
                .iter()
                .map(ProfileLoop::polygon_points)
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    // Ordinary nondegenerate one-loop lofts can be nonplanar, and retain
    // the general analytic homotopy skin. Only an actual correspondence
    // collapse or several loops calls for the planar limiting construction.
    let needed = polygons.iter().any(|p| {
        p.len() > 1
            || p.iter().any(|ring| {
                ring.iter()
                    .zip(ring.iter().cycle().skip(1))
                    .any(|(&a, &b)| key(a) == key(b))
            })
    });
    if !needed {
        return Ok(None);
    }
    orient_polygons(sections, &mut polygons)?;
    let mut skin = Skin::default();
    let topology = sections
        .iter()
        .zip(&polygons)
        .map(|(s, p)| skin.section(s, p))
        .collect::<Vec<_>>();
    let mut shell = Shell::new();
    skin.cap(&mut shell, &topology[0], true)?;
    for pair in topology.windows(2) {
        skin.connect(&mut shell, &pair[0], &pair[1])?;
    }
    skin.cap(
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
            shell.iter().flat_map(|f| f.edge_iter()).all(|e| e
                .front()
                .point()
                .distance(e.back().point())
                > 0.0)
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
