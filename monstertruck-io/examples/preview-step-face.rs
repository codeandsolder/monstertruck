//! Renders selected STEP faces to diagnostic preview images.

use anyhow::{Context, Result, anyhow, bail};
use clap::Parser;
use image::{Rgb, RgbImage};
use monstertruck_io::step::load::{
    Table,
    step_geometry::{
        Attributes, BoundedCurve, BoundingBox, Curve3D, EuclideanSpace, InnerSpace, Invertible,
        Origin, ParametricCurve3D, ParametricSurface3D, Point2, Point3, PolygonMesh, PolylineCurve,
        StepParameterCurve, Surface, TOLERANCE, Vector2, Vector3, Zero,
    },
};
use monstertruck_meshing::prelude::*;
use monstertruck_topology::compress::{CompressedEdgeUse, CompressedShell, CompressedTrimmedShell};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

type MeshShell = CompressedShell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>;

const BACKGROUND: Rgb<u8> = Rgb([30, 33, 38]);
const FILL: Rgb<u8> = Rgb([220, 188, 80]);
const WIRE: Rgb<u8> = Rgb([44, 132, 255]);
const BOUNDARY: Rgb<u8> = Rgb([255, 218, 0]);

#[derive(Parser, Debug)]
struct Args {
    /// Input STEP file.
    input_step_file: PathBuf,
    /// One-based face indices to render.
    #[arg(long, required = true)]
    face: Vec<usize>,
    /// One-based shell index to render.
    #[arg(long, default_value_t = 1)]
    shell: usize,
    /// Output directory.
    #[arg(long, default_value = "target/face-previews")]
    out: PathBuf,
    /// Preview image size in pixels.
    #[arg(long, default_value_t = 1200)]
    size: u32,
    /// Relative tolerance factor multiplied by shell bounding-box diameter.
    #[arg(long, default_value_t = 0.001)]
    tolerance_factor: f64,
    /// Print raw face trim ranges before tessellation.
    #[arg(long)]
    dump_trims: bool,
}

#[derive(Clone, Copy)]
struct Bounds2 {
    min: Vector2,
    max: Vector2,
}

impl Bounds2 {
    fn from_points(points: impl IntoIterator<Item = Vector2>) -> Option<Self> {
        points.into_iter().fold(None, |bounds, point| {
            Some(bounds.map_or(
                Self {
                    min: point,
                    max: point,
                },
                |bounds| Self {
                    min: Vector2::new(bounds.min.x.min(point.x), bounds.min.y.min(point.y)),
                    max: Vector2::new(bounds.max.x.max(point.x), bounds.max.y.max(point.y)),
                },
            ))
        })
    }

    #[expect(
        clippy::suboptimal_flops,
        reason = "diagnostic pixel projection keeps explicit multiply/add rounding; FMA changes pixel rounding at boundaries"
    )]
    fn project(self, point: Vector2, size: u32) -> (i32, i32) {
        let size = f64::from(size);
        let padding = 0.06 * size;
        let width = (self.max.x - self.min.x).max(TOLERANCE);
        let height = (self.max.y - self.min.y).max(TOLERANCE);
        let drawable = size - padding * 2.0;
        let scale = (drawable / width).min(drawable / height);
        let center = (self.min + self.max) * 0.5;
        let x = (point.x - center.x) * scale + size * 0.5;
        let y = size * 0.5 - (point.y - center.y) * scale;
        (rounded_i32(x), rounded_i32(y))
    }
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "value is rounded and clamped to the complete i32 range before the final audited float-to-integer cast"
)]
fn rounded_i32(value: f64) -> i32 {
    if value.is_nan() {
        return 0;
    }
    let rounded = value
        .round()
        .clamp(f64::from(i32::MIN), f64::from(i32::MAX));
    rounded as i32
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.size == 0 {
        bail!("preview size must be nonzero");
    }
    let step_file = fs::read_to_string(&args.input_step_file)
        .with_context(|| format!("failed to read {}", args.input_step_file.display()))?;
    let table = Table::from_step(&step_file).context("failed to parse STEP")?;
    let shell = load_shell(&table, args.shell)?;
    let tolerance = f64::max(
        shell_bounding_box(&shell).diameter() * args.tolerance_factor,
        TOLERANCE,
    );
    if args.dump_trims {
        args.face
            .iter()
            .try_for_each(|face| dump_trim_boundaries(&shell, *face, tolerance))?;
    }
    let meshed = shell.robust_triangulation(tolerance);
    fs::create_dir_all(&args.out)
        .with_context(|| format!("failed to create {}", args.out.display()))?;
    args.face
        .into_iter()
        .try_for_each(|face| render_face_pair(&meshed, face, &args.out, args.size))
}

fn load_shell(
    table: &Table,
    one_based_shell: usize,
) -> Result<CompressedTrimmedShell<Point3, Curve3D, Surface, StepParameterCurve>> {
    if one_based_shell == 0 {
        bail!("shell indices are one-based");
    }
    let mut shell_entries = table.shell.iter().collect::<Vec<_>>();
    shell_entries.sort_by_key(|(id, _)| *id);
    let (_, shell_holder) = shell_entries
        .get(one_based_shell - 1)
        .ok_or_else(|| anyhow!("shell {one_based_shell} not found"))?;
    table
        .to_compressed_trimmed_shell(shell_holder)
        .map_err(|error| anyhow!("failed to convert STEP shell: {error}"))
}

fn shell_bounding_box<C, S, T>(
    shell: &CompressedTrimmedShell<Point3, C, S, T>,
) -> BoundingBox<Point3>
where
    C: ParametricCurve3D + BoundedCurve,
    S: ParametricSurface3D, {
    let mut bounds: BoundingBox<Point3> = shell.vertices.iter().collect();
    for edge in &shell.edges {
        let (front, back) = edge.curve.range_tuple();
        for index in 0..=4 {
            let parameter = front + (back - front) * f64::from(index) / 4.0;
            bounds.push(edge.curve.evaluate(parameter));
        }
    }
    bounds
}

fn dump_trim_boundaries(
    shell: &CompressedTrimmedShell<Point3, Curve3D, Surface, StepParameterCurve>,
    one_based_face: usize,
    tolerance: f64,
) -> Result<()> {
    if one_based_face == 0 {
        bail!("face indices are one-based");
    }
    let face = shell
        .faces
        .get(one_based_face - 1)
        .ok_or_else(|| anyhow!("face {one_based_face} not found"))?;
    let orientation = face.orientation;
    let loop_count = face.boundaries.len();
    eprintln!("trim_dump face={one_based_face} orientation={orientation} loops={loop_count}");
    face.boundaries.iter().enumerate().for_each(|(loop_index, wire)| {
        let edge_uses = oriented_edge_uses(wire, face.orientation);
        edge_uses.into_iter().enumerate().for_each(
            |(edge_position, (edge_index, orientation, trim_curve))| {
                let points = trim_curve
                    .as_ref()
                    .map_or_default(|trim_curve| trim_curve.exact_trim_boundary_2d(tolerance));
                let bounds = Bounds2::from_points(
                    points.iter().map(|point| Vector2::new(point.x, point.y)),
                );
                let front = points
                    .first()
                    .map_or_else(|| "-".to_string(), format_point);
                let back = points
                    .last()
                    .map_or_else(|| "-".to_string(), format_point);
                let range = bounds.map_or_else(
                    || "-".to_string(),
                    |bounds| {
                        format!(
                            "({:.6},{:.6})..({:.6},{:.6})",
                            bounds.min.x,
                            bounds.min.y,
                            bounds.max.x,
                            bounds.max.y,
                        )
                    },
                );
                let point_count = points.len();
                eprintln!(
                    "trim_dump face={one_based_face} loop={loop_index} edge_pos={edge_position} edge={edge_index} orientation={orientation} points={point_count} front={front} back={back} range={range}",
                );
            },
        );
    });
    Ok(())
}

fn oriented_edge_uses<T>(
    wire: &[CompressedEdgeUse<T>],
    face_orientation: bool,
) -> Vec<(usize, bool, Option<&T>)> {
    if face_orientation {
        wire.iter()
            .map(|edge_use| {
                (
                    edge_use.index,
                    edge_use.orientation,
                    edge_use.trim_curve.as_ref(),
                )
            })
            .collect()
    } else {
        wire.iter()
            .rev()
            .map(|edge_use| {
                (
                    edge_use.index,
                    !edge_use.orientation,
                    edge_use.trim_curve.as_ref(),
                )
            })
            .collect()
    }
}

fn format_point(point: &Point2) -> String { format!("({:.6},{:.6})", point.x, point.y) }

fn render_face_pair(
    shell: &MeshShell,
    one_based_face: usize,
    output_dir: &Path,
    size: u32,
) -> Result<()> {
    if one_based_face == 0 {
        bail!("face indices are one-based");
    }
    let face = shell
        .faces
        .get(one_based_face - 1)
        .ok_or_else(|| anyhow!("face {one_based_face} not found"))?;
    let mesh = face
        .surface
        .as_ref()
        .ok_or_else(|| anyhow!("face {one_based_face} did not mesh"))?;
    let mesh = if face.orientation {
        mesh.clone()
    } else {
        mesh.inverse()
    };
    let prefix = format!("shell-001-face-{one_based_face:03}");
    render_mesh_uv(&mesh, &output_dir.join(format!("{prefix}-uv.png")), size)?;
    render_mesh_3d(&mesh, &output_dir.join(format!("{prefix}-3d.png")), size)
}

fn render_mesh_uv(mesh: &PolygonMesh, path: &Path, size: u32) -> Result<()> {
    let uv_coords = mesh.uv_coords();
    if uv_coords.is_empty() {
        bail!("mesh has no UV coordinates");
    }
    let points = uv_coords
        .iter()
        .map(|point| Vector2::new(point.x, point.y))
        .collect::<Vec<_>>();
    let bounds = Bounds2::from_points(points.iter().copied())
        .ok_or_else(|| anyhow!("mesh has no UV points"))?;
    let project = |index: usize| bounds.project(points[index], size);
    let triangles = mesh
        .tri_faces()
        .iter()
        .map(|triangle| {
            [
                project(triangle[0].uv.unwrap_or(triangle[0].pos)),
                project(triangle[1].uv.unwrap_or(triangle[1].pos)),
                project(triangle[2].uv.unwrap_or(triangle[2].pos)),
            ]
        })
        .collect::<Vec<_>>();
    write_preview(path, size, &triangles)
}

fn render_mesh_3d(mesh: &PolygonMesh, path: &Path, size: u32) -> Result<()> {
    let positions = mesh.positions();
    let (origin, axis_x, axis_y) = view_frame(mesh)?;
    let points = positions
        .iter()
        .map(|point| {
            let offset = *point - origin;
            Vector2::new(offset.dot(axis_x), offset.dot(axis_y))
        })
        .collect::<Vec<_>>();
    let bounds = Bounds2::from_points(points.iter().copied())
        .ok_or_else(|| anyhow!("mesh has no positions"))?;
    let project = |index: usize| bounds.project(points[index], size);
    let triangles = mesh
        .tri_faces()
        .iter()
        .map(|triangle| {
            [
                project(triangle[0].pos),
                project(triangle[1].pos),
                project(triangle[2].pos),
            ]
        })
        .collect::<Vec<_>>();
    write_preview(path, size, &triangles)
}

fn view_frame(mesh: &PolygonMesh) -> Result<(Point3, Vector3, Vector3)> {
    let positions = mesh.positions();
    if positions.is_empty() {
        bail!("cannot derive view frame from an empty mesh");
    }
    let origin = positions
        .iter()
        .copied()
        .fold(Point3::origin(), |sum, point| {
            Point3::from_vec(sum.to_vec() + point.to_vec())
        });
    let position_count =
        u32::try_from(positions.len()).context("preview mesh has more than u32::MAX positions")?;
    let origin = Point3::from_vec(origin.to_vec() / f64::from(position_count));
    let normal = mesh
        .normals()
        .iter()
        .copied()
        .fold(Vector3::zero(), |sum, normal| sum + normal);
    let normal = if normal.so_small() {
        mesh.tri_faces()
            .iter()
            .find_map(|triangle| {
                let a = positions[triangle[0].pos];
                let b = positions[triangle[1].pos];
                let c = positions[triangle[2].pos];
                let normal = (b - a).cross(c - a);
                (!normal.so_small()).then_some(normal.normalize())
            })
            .ok_or_else(|| anyhow!("cannot derive view normal"))?
    } else {
        normal.normalize()
    };
    let seed = if normal.x.abs() < 0.8 {
        Vector3::unit_x()
    } else {
        Vector3::unit_y()
    };
    let axis_x = normal.cross(seed).normalize();
    let axis_y = normal.cross(axis_x).normalize();
    Ok((origin, axis_x, axis_y))
}

fn write_preview(path: &Path, size: u32, triangles: &[[(i32, i32); 3]]) -> Result<()> {
    let mut image = RgbImage::from_pixel(size, size, BACKGROUND);
    for &triangle in triangles {
        fill_triangle(&mut image, triangle);
    }
    let boundary_edges = boundary_edges(triangles);
    for &triangle in triangles {
        draw_line(&mut image, triangle[0], triangle[1], WIRE);
        draw_line(&mut image, triangle[1], triangle[2], WIRE);
        draw_line(&mut image, triangle[2], triangle[0], WIRE);
    }
    for (front, back) in boundary_edges {
        draw_thick_line(&mut image, front, back, BOUNDARY);
    }
    image
        .save(path)
        .with_context(|| format!("failed to write {}", path.display()))
}

fn boundary_edges(triangles: &[[(i32, i32); 3]]) -> Vec<((i32, i32), (i32, i32))> {
    let mut counts = BTreeMap::<((i32, i32), (i32, i32)), usize>::new();
    for triangle in triangles {
        for (front, back) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            let key = if front <= back {
                (front, back)
            } else {
                (back, front)
            };
            *counts.entry(key).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .filter_map(|(edge, count)| (count == 1).then_some(edge))
        .collect()
}

fn fill_triangle(image: &mut RgbImage, triangle: [(i32, i32); 3]) {
    let Ok(width) = i32::try_from(image.width()) else {
        return;
    };
    let Ok(height) = i32::try_from(image.height()) else {
        return;
    };
    if width == 0 || height == 0 {
        return;
    }

    let [a, b, c] = triangle;
    let min_x = a.0.min(b.0).min(c.0).clamp(0, width - 1);
    let max_x = a.0.max(b.0).max(c.0).clamp(0, width - 1);
    let min_y = a.1.min(b.1).min(c.1).clamp(0, height - 1);
    let max_y = a.1.max(b.1).max(c.1).clamp(0, height - 1);
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            if point_in_triangle((x, y), triangle) {
                image.put_pixel(x.cast_unsigned(), y.cast_unsigned(), FILL);
            }
        }
    }
}

fn point_in_triangle(point: (i32, i32), triangle: [(i32, i32); 3]) -> bool {
    let edge = |a: (i32, i32), b: (i32, i32), p: (i32, i32)| {
        let px = i64::from(p.0);
        let py = i64::from(p.1);
        let ax = i64::from(a.0);
        let ay = i64::from(a.1);
        let bx = i64::from(b.0);
        let by = i64::from(b.1);
        (px - ax) * (by - ay) - (py - ay) * (bx - ax)
    };
    let [a, b, c] = triangle;
    let w0 = edge(b, c, point);
    let w1 = edge(c, a, point);
    let w2 = edge(a, b, point);
    (w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0)
}

fn draw_thick_line(image: &mut RgbImage, front: (i32, i32), back: (i32, i32), color: Rgb<u8>) {
    for dy in -1..=1 {
        for dx in -1..=1 {
            draw_line(
                image,
                (front.0.saturating_add(dx), front.1.saturating_add(dy)),
                (back.0.saturating_add(dx), back.1.saturating_add(dy)),
                color,
            );
        }
    }
}

fn draw_line(image: &mut RgbImage, front: (i32, i32), back: (i32, i32), color: Rgb<u8>) {
    let width = i64::from(image.width());
    let height = i64::from(image.height());
    let (mut x0, mut y0) = (i64::from(front.0), i64::from(front.1));
    let (x1, y1) = (i64::from(back.0), i64::from(back.1));
    let dx = (x1 - x0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let dy = -(y1 - y0).abs();
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut error = dx + dy;
    loop {
        if x0 >= 0
            && y0 >= 0
            && x0 < width
            && y0 < height
            && let (Ok(x), Ok(y)) = (u32::try_from(x0), u32::try_from(y0))
        {
            image.put_pixel(x, y, color);
        }
        if x0 == x1 && y0 == y1 {
            break;
        }
        let error2 = 2 * error;
        if error2 >= dy {
            error += dy;
            x0 += sx;
        }
        if error2 <= dx {
            error += dx;
            y0 += sy;
        }
    }
}
