//! Compare standard and robust tessellation of a compressed native shell.
//! Reports missing faces by carrier and topology after component-local welding.

use anyhow::{Context, Result};
use monstertruck_geometry::prelude::Invertible;
use monstertruck_mesh::{KERNEL_EPSILON, PolygonMesh};
use monstertruck_meshing::{
    analyzers::Topology,
    filters::OptimizingFilter,
    tessellation::{MeshableShape, MeshedShape, RobustMeshableShape},
};
use monstertruck_modeling::{
    Point3,
    geometry::{Curve, Surface},
};
use monstertruck_topology::compress::CompressedShell;
use serde_json::{Value, json};
use std::{collections::BTreeMap, fs::File, io::BufReader, path::PathBuf};

type CShell = CompressedShell<Point3, Curve, Surface>;

const fn kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane(_) => "plane",
        Surface::BsplineSurface(_) => "bspline",
        Surface::NurbsSurface(_) => "nurbs",
        Surface::Revolution(_) => "revolution",
        Surface::Tspline(_) => "tspline",
        Surface::Sphere(_) => "sphere",
        Surface::Torus(_) => "torus",
    }
}

fn stats(mut polygon: PolygonMesh) -> Value {
    polygon
        .put_together_same_attrs(KERNEL_EPSILON * 10.0)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    json!({"condition":format!("{:?}",polygon.shell_condition()),
        "polygons":polygon.faces().len(),"vertices":polygon.positions().len()})
}

// Face components follow the original native edge identities, not
// coincident coordinates. This deliberately does not weld between them.
fn components(shell: &CShell) -> Vec<Vec<usize>> {
    let mut owners = vec![None::<usize>; shell.edges.len()];
    let mut adjacency = vec![Vec::new(); shell.faces.len()];
    for (face_index, face) in shell.faces.iter().enumerate() {
        for edge in face.boundaries.iter().flatten() {
            match owners[edge.index] {
                Some(other) => {
                    adjacency[face_index].push(other);
                    adjacency[other].push(face_index);
                }
                None => owners[edge.index] = Some(face_index),
            }
        }
    }
    let mut seen = vec![false; shell.faces.len()];
    let mut groups = Vec::new();
    for start in 0..seen.len() {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut pending = vec![start];
        let mut group = Vec::new();
        while let Some(face) = pending.pop() {
            group.push(face);
            for &neighbor in &adjacency[face] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    pending.push(neighbor);
                }
            }
        }
        group.sort_unstable();
        groups.push(group);
    }
    groups
}

fn inspect(shell: &CShell, tolerance: f64, robust: bool, groups: &[Vec<usize>]) -> Value {
    let mesh = if robust {
        shell.robust_triangulation(tolerance)
    } else {
        shell.triangulation(tolerance)
    };
    let mut missing = BTreeMap::<&str, Vec<usize>>::new();
    let mut empty = BTreeMap::<&str, Vec<usize>>::new();
    for (i, (face, tessellated)) in shell.faces.iter().zip(&mesh.faces).enumerate() {
        match &tessellated.surface {
            None => missing.entry(kind(&face.surface)).or_default().push(i),
            Some(polygon) if polygon.faces().is_empty() => {
                empty.entry(kind(&face.surface)).or_default().push(i)
            }
            Some(_) => {}
        }
    }
    let components = groups
        .iter()
        .map(|group| {
            let mut polygon = PolygonMesh::default();
            for &index in group {
                let face = &mesh.faces[index];
                if let Some(surface) = &face.surface {
                    if face.orientation {
                        polygon.merge(surface.clone());
                    } else {
                        polygon.merge(surface.inverse());
                    }
                }
            }
            json!({"native_faces": group.len(), "polygon":stats(polygon)})
        })
        .collect::<Vec<_>>();
    json!({"native_faces":shell.faces.len(),"missing":missing,"empty":empty,
        "polygon":stats(mesh.to_polygon()),"components":components})
}

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let path = PathBuf::from(args.next().context("expected compressed shell JSON")?);
    let tolerance: f64 = args
        .next()
        .context("expected tessellation tolerance")?
        .to_str()
        .context("non-UTF8 tolerance")?
        .parse()?;
    anyhow::ensure!(
        args.next().is_none() && tolerance.is_finite() && tolerance > 0.0,
        "usage: inspect-tessellation SHELL_JSON TOLERANCE"
    );
    let shell: CShell = serde_json::from_reader(BufReader::new(File::open(path)?))?;
    let groups = components(&shell);
    let standard = inspect(&shell, tolerance, false, &groups);
    let robust = inspect(&shell, tolerance, true, &groups);
    println!("{}", json!({"standard":standard,"robust":robust}));
    Ok(())
}
