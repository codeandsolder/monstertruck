//! Profiles STEP tessellation at an explicit tolerance factor.

use anyhow::{Context, Result};
use clap::Parser;
use monstertruck_io::step::load::{
    Table,
    convert::StepCompressedTrimmedShell,
    step_geometry::{
        BoundedCurve, Curve3D, ElementarySurface, InnerSpace, MetricSpace, ParameterBoundary2D,
        ParameterDivision1D, Point3, PolygonMesh, PolylineCurve, SearchNearestParameter,
        SearchParameter, Surface, SweepSurface,
    },
};
use monstertruck_meshing::prelude::*;
use monstertruck_topology::compress::CompressedShell;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
struct Args {
    /// Input STEP file.
    input_step_file: PathBuf,
    /// Relative tolerance factor multiplied by shell bbox diameter.
    #[arg(long, default_value_t = 0.001)]
    tolerance_factor: f64,
    /// Optional shell index filter.
    #[arg(long)]
    shell: Option<usize>,
    /// Optional face indices to inspect before meshing.
    #[arg(long, value_delimiter = ',')]
    inspect_faces: Vec<usize>,
    /// Tessellation mode.
    #[arg(long, default_value = "robust")]
    mode: TessellationMode,
    /// Report shared-edge vs face-boundary drift after meshing.
    #[arg(long)]
    check_boundaries: bool,
}

#[derive(clap::ValueEnum, Clone, Copy, Debug)]
enum TessellationMode {
    Plain,
    Robust,
}

const fn curve_kind(curve: &Curve3D) -> &'static str {
    match curve {
        Curve3D::Line(_) => "Line",
        Curve3D::Polyline(_) => "Polyline",
        Curve3D::Conic(_) => "Conic",
        Curve3D::BsplineCurve(_) => "BsplineCurve",
        Curve3D::ParameterCurve(_) => "ParameterCurve",
        Curve3D::SurfaceCurve(_) => "SurfaceCurve",
        Curve3D::IntersectionCurve(_) => "IntersectionCurve",
        Curve3D::NurbsCurve(_) => "NurbsCurve",
    }
}

const fn surface_kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Elementary(ElementarySurface::Plane(_)) => "Plane",
        Surface::Elementary(ElementarySurface::Sphere(_)) => "Sphere",
        Surface::Elementary(ElementarySurface::CylindricalSurface(_)) => "CylindricalSurface",
        Surface::Elementary(ElementarySurface::ToroidalSurface(_)) => "ToroidalSurface",
        Surface::Elementary(ElementarySurface::ConicalSurface(_)) => "ConicalSurface",
        Surface::Sweep(SweepSurface::ExtrusionSurface(_)) => "ExtrusionSurface",
        Surface::Sweep(SweepSurface::RevolutionSurface(_)) => "RevolutionSurface",
        Surface::Bspline(_) => "BsplineSurface",
        Surface::Nurbs(_) => "NurbsSurface",
    }
}

fn point_segment_distance(point: Point3, front: Point3, back: Point3) -> f64 {
    let segment = back - front;
    let denom = segment.dot(segment);
    if denom.so_small() {
        point.distance(front)
    } else {
        let t = ((point - front).dot(segment) / denom).clamp(0.0, 1.0);
        let nearest = front + segment * t;
        point.distance(nearest)
    }
}

fn mesh_boundary_segments(mesh: &PolygonMesh) -> Vec<(Point3, Point3)> {
    let positions = mesh.positions();
    let mut edge_counts = std::collections::BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
    for tri in mesh.tri_faces() {
        let indices = [tri[0].pos, tri[1].pos, tri[2].pos];
        for i in 0..3 {
            let front = indices[i];
            let back = indices[(i + 1) % 3];
            let key = if front < back {
                (front, back)
            } else {
                (back, front)
            };
            edge_counts.entry(key).or_default().push((front, back));
        }
    }
    edge_counts
        .into_iter()
        .filter(|(_, occurrences)| occurrences.len() == 1)
        .map(|(_, occurrences)| occurrences[0])
        .map(|(front, back)| (positions[front], positions[back]))
        .collect()
}

fn report_boundary_drift(
    shell_idx: usize,
    shell: &CompressedShell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
) {
    let face_segments = shell
        .faces
        .iter()
        .map(|face| {
            face.surface
                .as_ref()
                .map_or_else(Vec::new, mesh_boundary_segments)
        })
        .collect::<Vec<_>>();
    for (edge_idx, edge) in shell.edges.iter().enumerate() {
        let incident_faces = shell
            .faces
            .iter()
            .enumerate()
            .filter(|(_, face)| {
                face.boundaries
                    .iter()
                    .flatten()
                    .any(|boundary_edge| boundary_edge.index == edge_idx)
            })
            .map(|(face_idx, _)| face_idx)
            .collect::<Vec<_>>();
        let relevant_segments = shell
            .faces
            .iter()
            .enumerate()
            .filter(|(face_idx, _)| incident_faces.contains(face_idx))
            .flat_map(|(face_idx, _)| face_segments[face_idx].iter().copied())
            .collect::<Vec<_>>();
        if relevant_segments.is_empty() {
            continue;
        }
        let max_drift = edge
            .curve
            .iter()
            .copied()
            .map(|point| {
                relevant_segments
                    .iter()
                    .map(|(front, back)| point_segment_distance(point, *front, *back))
                    .min_by(f64::total_cmp)
                    .unwrap_or(f64::INFINITY)
            })
            .max_by(f64::total_cmp)
            .unwrap_or(0.0);
        if max_drift > 1.0e-6 {
            eprintln!(
                "boundary-drift shell={shell_idx} edge={edge_idx} edge_points={} incident_faces={incident_faces:?} segment_count={} max_drift={max_drift:.9}",
                edge.curve.len(),
                relevant_segments.len(),
            );
        }
    }
}

fn inspect_shell_faces(
    shell_idx: usize,
    shell: &StepCompressedTrimmedShell,
    inspect_faces: &[usize],
    tolerance: f64,
) {
    for &face_idx in inspect_faces {
        let Some(face) = shell.faces.get(face_idx) else {
            continue;
        };
        eprintln!(
            "inspect shell={shell_idx} face={face_idx} surface={} loops={}",
            surface_kind(&face.surface),
            face.boundaries.len(),
        );
        for (wire_idx, wire) in face.boundaries.iter().enumerate() {
            let kinds = wire
                .iter()
                .filter_map(|edge_idx| shell.edges.get(edge_idx.index))
                .map(|edge| curve_kind(&edge.curve))
                .collect::<Vec<_>>();
            eprintln!(
                "inspect shell={shell_idx} face={face_idx} wire={wire_idx} edge_kinds={kinds:?}",
            );
            for (edge_pos, edge_idx) in wire.iter().enumerate() {
                let Some(edge) = shell.edges.get(edge_idx.index) else {
                    continue;
                };
                let boundary = edge.curve.parameter_boundary_2d(&face.surface, tolerance);
                let boundary_points = boundary.as_ref().map_or(0, Vec::len);
                let boundary_front = boundary
                    .as_ref()
                    .and_then(|points| points.first().copied())
                    .map_or_else(|| "-".into(), |uv| format!("({:.6},{:.6})", uv.x, uv.y));
                let boundary_back = boundary
                    .as_ref()
                    .and_then(|points| points.last().copied())
                    .map_or_else(|| "-".into(), |uv| format!("({:.6},{:.6})", uv.x, uv.y));
                let edge_curve_points = edge
                    .curve
                    .parameter_division(edge.curve.range_tuple(), tolerance)
                    .1
                    .len();
                match &edge.curve {
                    Curve3D::SurfaceCurve(curve) => {
                        let edge_index = edge_idx.index;
                        let exact = boundary.is_some();
                        let assoc = curve.associated_geometry().len();
                        let face_surface = surface_kind(&face.surface);
                        eprintln!(
                            "inspect shell={shell_idx} face={face_idx} wire={wire_idx} edge={edge_pos} edge_index={edge_index} exact={exact} boundary_points={boundary_points} edge_curve_points={edge_curve_points} boundary_front={boundary_front} boundary_back={boundary_back} assoc={assoc} face_surface={face_surface}",
                        );
                    }
                    Curve3D::IntersectionCurve(curve) => {
                        let edge_index = edge_idx.index;
                        let exact = boundary.is_some();
                        let surface0 = surface_kind(curve.surface0().as_ref());
                        let surface1 = surface_kind(curve.surface1().as_ref());
                        let face_surface = surface_kind(&face.surface);
                        eprintln!(
                            "inspect shell={shell_idx} face={face_idx} wire={wire_idx} edge={edge_pos} edge_index={edge_index} exact={exact} boundary_points={boundary_points} edge_curve_points={edge_curve_points} boundary_front={boundary_front} boundary_back={boundary_back} s0={surface0} s1={surface1} face_surface={face_surface}",
                        );
                    }
                    _ => {
                        let edge_index = edge_idx.index;
                        let orientation = edge_idx.orientation;
                        let exact = boundary.is_some();
                        eprintln!(
                            "inspect shell={shell_idx} face={face_idx} wire={wire_idx} edge={edge_pos} edge_index={edge_index} orientation={orientation} exact={exact} boundary_points={boundary_points} edge_curve_points={edge_curve_points} boundary_front={boundary_front} boundary_back={boundary_back}",
                        );
                    }
                }
                if boundary.is_none() {
                    for (sample_idx, point) in edge
                        .curve
                        .parameter_division(edge.curve.range_tuple(), tolerance)
                        .1
                        .into_iter()
                        .enumerate()
                    {
                        let direct = face
                            .surface
                            .search_parameter(point, None, 100)
                            .map_or_else(|| "-".into(), |(u, v)| format!("({u:.6},{v:.6})"));
                        let nearest = face
                            .surface
                            .search_nearest_parameter(point, None, 100)
                            .map_or_else(|| "-".into(), |(u, v)| format!("({u:.6},{v:.6})"));
                        eprintln!(
                            "inspect shell={shell_idx} face={face_idx} wire={wire_idx} edge={edge_pos} sample={sample_idx} point=({:.6},{:.6},{:.6}) direct={direct} nearest={nearest}",
                            point.x, point.y, point.z,
                        );
                    }
                }
            }
        }
    }
}

fn main() -> Result<()> {
    let Args {
        input_step_file,
        tolerance_factor,
        shell,
        inspect_faces,
        mode,
        check_boundaries,
    } = Args::parse();
    let step_file = std::fs::read_to_string(&input_step_file)
        .with_context(|| format!("failed to read {}", input_step_file.display()))?;
    let table = Table::from_step(&step_file).context("failed to parse STEP")?;
    let assy = table.step_assy().context("failed to build STEP assembly")?;
    let mut shells = assy
        .top_nodes()
        .flat_map(|top| assy.paths_iter(top.index()))
        .flat_map(|path| {
            path.terminal_node()
                .into_iter()
                .flat_map(|node| node.shape().iter())
                .filter_map(|idx| {
                    if let Some(step_solid) = table.manifold_solid_brep.get(idx) {
                        table
                            .to_compressed_trimmed_solid(step_solid)
                            .ok()
                            .map(|solid| solid.boundaries)
                    } else if let Some(step_shells) = table.shell_based_surface_model.get(idx) {
                        table.to_compressed_trimmed_shells(step_shells).ok()
                    } else {
                        None
                    }
                })
                .flatten()
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    if shells.is_empty() {
        shells = table
            .shell
            .values()
            .filter_map(|shell| table.to_compressed_trimmed_shell(shell).ok())
            .collect();
    }
    shells.sort_by_key(|shell| shell.faces.len());

    for (shell_idx, shell_data) in shells.into_iter().enumerate() {
        if shell.is_some_and(|target| shell_idx != target) {
            continue;
        }
        let diameter = shell_data.bounding_box().diameter();
        let tolerance = shell_data.relative_tolerance(tolerance_factor);
        eprintln!(
            "profile shell={shell_idx} faces={} diameter={diameter:.6} tolerance={tolerance:.9}",
            shell_data.faces.len(),
        );
        inspect_shell_faces(shell_idx, &shell_data, &inspect_faces, tolerance);

        let start = Instant::now();
        let meshed = match mode {
            TessellationMode::Plain => shell_data.triangulation(tolerance),
            TessellationMode::Robust => shell_data.robust_triangulation(tolerance),
        };
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        let meshed_faces = meshed
            .faces
            .iter()
            .filter(|face| face.surface.is_some())
            .count();
        eprintln!(
            "profile shell={shell_idx} mode={mode:?} elapsed_ms={elapsed_ms:.3} meshed_faces={meshed_faces}/{}",
            meshed.faces.len(),
        );
        if check_boundaries {
            report_boundary_drift(shell_idx, &meshed);
        }
    }
    Ok(())
}
