//! Orthographic material-silhouette comparisons on independent OCCT meshes.
use super::compare::Thresholds;
use anyhow::{Context, Result, ensure};
use clap::ValueEnum;
use image::{GrayImage, Luma, imageops};
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

type Point = [f64; 3];
type Triangle = [Point; 3];
const SIZE: u32 = 320;
const SUPERSAMPLE: u32 = 2;
const MARGIN: f64 = 0.06;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
#[value(rename_all = "snake_case")]
pub enum View {
    Px,
    Nx,
    Py,
    Ny,
    Pz,
    Nz,
    IsoPpp,
    IsoNpp,
}

impl View {
    pub const ALL: [Self; 8] = [
        Self::Px,
        Self::Nx,
        Self::Py,
        Self::Ny,
        Self::Pz,
        Self::Nz,
        Self::IsoPpp,
        Self::IsoNpp,
    ];
    fn direction(self) -> Point {
        match self {
            Self::Px => [1.0, 0.0, 0.0],
            Self::Nx => [-1.0, 0.0, 0.0],
            Self::Py => [0.0, 1.0, 0.0],
            Self::Ny => [0.0, -1.0, 0.0],
            Self::Pz => [0.0, 0.0, 1.0],
            Self::Nz => [0.0, 0.0, -1.0],
            Self::IsoPpp => [1.0, 1.0, 1.0],
            Self::IsoNpp => [-1.0, 1.0, 1.0],
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::Px => "px",
            Self::Nx => "nx",
            Self::Py => "py",
            Self::Ny => "ny",
            Self::Pz => "pz",
            Self::Nz => "nz",
            Self::IsoPpp => "iso_ppp",
            Self::IsoNpp => "iso_npp",
        }
    }
}

fn dot(a: Point, b: Point) -> f64 { a.into_iter().zip(b).map(|(x, y)| x * y).sum() }
fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: Point, b: Point) -> Point { std::array::from_fn(|i| a[i] - b[i]) }
fn unit(p: Point) -> Point {
    let norm = dot(p, p).sqrt();
    p.map(|x| x / norm)
}

fn triangles(coordinates: &[f64], shift: Point) -> Result<Vec<Triangle>> {
    let (chunks, remainder) = coordinates.as_chunks::<9>();
    ensure!(remainder.is_empty(), "incomplete OCCT triangle");
    chunks
        .iter()
        .map(|coords| {
            let triangle =
                std::array::from_fn(|i| std::array::from_fn(|j| coords[3 * i + j] + shift[j]));
            ensure!(
                triangle.iter().flatten().all(|x: &f64| x.is_finite()),
                "nonfinite triangle"
            );
            Ok(triangle)
        })
        .collect()
}

#[derive(Debug, Clone, Copy)]
struct Frame {
    u: Point,
    v: Point,
    d: Point,
    xmin: f64,
    ymax: f64,
    sx: f64,
    sy: f64,
}

impl Frame {
    fn new(reference: &[Triangle], view: View) -> Result<Self> {
        ensure!(!reference.is_empty(), "reference mesh is empty");
        let d = unit(view.direction());
        let up = if d[2].abs() > 0.9 {
            [0.0, 1.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        let u = unit(cross(up, d));
        let v = unit(cross(d, u));
        let (mut xmin, mut xmax, mut ymin, mut ymax) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for &point in reference.iter().flatten() {
            let x = dot(point, u);
            let y = dot(point, v);
            xmin = xmin.min(x);
            xmax = xmax.max(x);
            ymin = ymin.min(y);
            ymax = ymax.max(y);
        }
        let width = (xmax - xmin).max(1e-9);
        let height = (ymax - ymin).max(1e-9);
        xmin -= width * MARGIN;
        xmax += width * MARGIN;
        ymin -= height * MARGIN;
        ymax += height * MARGIN;
        let span = f64::from(SIZE * SUPERSAMPLE - 1);
        let frame = Self {
            u,
            v,
            d,
            xmin,
            ymax,
            sx: span / (xmax - xmin).max(1e-12),
            sy: span / (ymax - ymin).max(1e-12),
        };
        ensure!(
            [frame.xmin, frame.ymax, frame.sx, frame.sy]
                .into_iter()
                .all(f64::is_finite),
            "nonfinite projection"
        );
        Ok(frame)
    }
    fn project(self, p: Point) -> Result<[i32; 2]> {
        let x = ((dot(p, self.u) - self.xmin) * self.sx).round_ties_even();
        let y = ((self.ymax - dot(p, self.v)) * self.sy).round_ties_even();
        Ok([
            x.to_i32()
                .context("projected x coordinate exceeds raster range")?,
            y.to_i32()
                .context("projected y coordinate exceeds raster range")?,
        ])
    }
}

struct Projected {
    points: [[i32; 2]; 3],
    depth: f64,
    value: u8,
}
impl Projected {
    fn new(triangle: Triangle, frame: Frame) -> Result<Self> {
        let mut normal = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
        let length = dot(normal, normal).sqrt();
        if length > 1e-18 {
            normal = normal.map(|v| v / length);
        }
        let shade = (0.28 + 0.67 * dot(normal, frame.d).abs()).clamp(0.0, 1.0);
        let depth = triangle.into_iter().map(|p| dot(p, frame.d)).sum::<f64>() / 3.0;
        ensure!(depth.is_finite(), "nonfinite projected depth");
        Ok(Self {
            points: [
                frame.project(triangle[0])?,
                frame.project(triangle[1])?,
                frame.project(triangle[2])?,
            ],
            value: (235.0 - 190.0 * shade)
                .round_ties_even()
                .to_u8()
                .context("nonfinite shade")?,
            depth,
        })
    }
}

// Inclusive, clipped triangle scanlines at supersampled integer coordinates.
// Convex triangles need only their minimum/maximum edge intersections.
// The opposing half-pixel tie rules avoid double-filling a half-pixel edge.
fn fill_triangle(image: &mut GrayImage, triangle: &Projected) -> Result<()> {
    let p = triangle.points;
    let width = i32::try_from(image.width())?;
    let height = i32::try_from(image.height())?;
    let ymin = p
        .iter()
        .map(|v| v[1])
        .min()
        .context("triangle has no points")?
        .max(0);
    let ymax = p
        .iter()
        .map(|v| v[1])
        .max()
        .context("triangle has no points")?
        .min(height - 1);
    let stride = usize::try_from(width)?;
    for y in ymin..=ymax {
        let (mut left, mut right) = (f64::INFINITY, f64::NEG_INFINITY);
        for (a, b) in [(p[0], p[1]), (p[1], p[2]), (p[2], p[0])] {
            if !(a[1].min(b[1])..=a[1].max(b[1])).contains(&y) {
                continue;
            }
            if a[1] == b[1] {
                left = left.min(f64::from(a[0].min(b[0])));
                right = right.max(f64::from(a[0].max(b[0])));
            } else {
                let x = f64::from(a[0])
                    + (f64::from(b[0]) - f64::from(a[0])) * (f64::from(y) - f64::from(a[1]))
                        / (f64::from(b[1]) - f64::from(a[1]));
                left = left.min(x);
                right = right.max(x);
            }
        }
        let first = (left + 0.5).floor().max(0.0);
        let last = (right - 0.5).ceil().min(f64::from(width - 1));
        if first <= last {
            let row = usize::try_from(y)? * stride;
            let from = first.to_usize().context("invalid left pixel")?;
            let to = last.to_usize().context("invalid right pixel")?;
            image
                .as_mut()
                .get_mut(row + from..=row + to)
                .context("invalid raster span")?
                .fill(triangle.value);
        }
    }
    Ok(())
}

fn render(triangles: &[Triangle], frame: Frame) -> Result<GrayImage> {
    let mut projected: Vec<_> = triangles
        .iter()
        .map(|&triangle| Projected::new(triangle, frame))
        .collect::<Result<_>>()?;
    // Equal depths retain mesh order, independent of sort implementation.
    projected.sort_by(|a, b| a.depth.total_cmp(&b.depth));
    let extent = SIZE * SUPERSAMPLE;
    let mut image = GrayImage::from_pixel(extent, extent, Luma([255]));
    for triangle in &projected {
        fill_triangle(&mut image, triangle)?;
    }
    Ok(imageops::resize(
        &image,
        SIZE,
        SIZE,
        imageops::FilterType::Lanczos3,
    ))
}

#[derive(Debug, Serialize)]
pub struct RenderComparison {
    view: View,
    mean_abs: f64,
    changed_fraction: f64,
    silhouette_xor: f64,
    reference: PathBuf,
    candidate: PathBuf,
    diff: PathBuf,
    pub pass: bool,
}

fn compare_pixels(reference: &GrayImage, candidate: &GrayImage) -> Result<(GrayImage, [f64; 3])> {
    ensure!(
        reference.dimensions() == candidate.dimensions() && !reference.is_empty(),
        "image dimension mismatch"
    );
    let (mut sum, mut changed, mut silhouette) = (0_u64, 0_u64, 0_u64);
    let mut differences = Vec::with_capacity(reference.len());
    for (&a, &b) in reference.as_raw().iter().zip(candidate.as_raw()) {
        let delta = a.abs_diff(b);
        sum += u64::from(delta);
        changed += u64::from(delta > 8);
        silhouette += u64::from((a < 250) != (b < 250));
        differences.push(delta);
    }
    let pixels = f64::from(reference.width()) * f64::from(reference.height());
    let diff = GrayImage::from_raw(reference.width(), reference.height(), differences)
        .context("invalid difference buffer")?;
    Ok((
        diff,
        [
            sum.to_f64().context("pixel sum overflow")? / pixels / 255.0,
            changed.to_f64().context("changed pixel count overflow")? / pixels,
            silhouette
                .to_f64()
                .context("silhouette pixel count overflow")?
                / pixels,
        ],
    ))
}

pub fn compare(
    reference: &[f64],
    candidate: &[f64],
    shift: Point,
    thresholds: &Thresholds,
    views: &[View],
    directory: &Path,
) -> Result<Vec<RenderComparison>> {
    let reference = triangles(reference, [0.0; 3])?;
    let candidate = triangles(candidate, shift)?;
    let limits = [
        thresholds.render_mean_abs.unwrap_or(1.0),
        thresholds.render_changed_fraction.unwrap_or(1.0),
        thresholds.silhouette_xor.unwrap_or(1.0),
    ];
    ensure!(
        limits.into_iter().all(|x| x.is_finite() && x >= 0.0),
        "invalid render threshold"
    );
    let views = if views.is_empty() {
        &View::ALL[..]
    } else {
        views
    };
    ensure!(
        views
            .iter()
            .enumerate()
            .all(|(i, v)| !views[..i].contains(v)),
        "duplicate render view"
    );
    std::fs::create_dir_all(directory)?;
    let mut results = Vec::with_capacity(views.len());
    // Contact sheet row order matches the reports, with 3 columns: reference,
    // candidate, absolute difference.
    let thumb = 220_u32;
    let mut sheet =
        GrayImage::from_pixel(thumb * 3, thumb * u32::try_from(views.len())?, Luma([255]));
    for (row, &view) in views.iter().enumerate() {
        let frame = Frame::new(&reference, view)?;
        let ref_img = render(&reference, frame)?;
        let cand_img = render(&candidate, frame)?;
        let (diff, metrics) = compare_pixels(&ref_img, &cand_img)?;
        let name = view.name();
        let ref_path = directory.join(format!("{name}.reference.png"));
        let cand_path = directory.join(format!("{name}.candidate.png"));
        let diff_path = directory.join(format!("{name}.diff.png"));
        ref_img.save(&ref_path)?;
        cand_img.save(&cand_path)?;
        diff.save(&diff_path)?;
        for (column, im) in [&ref_img, &cand_img, &diff].into_iter().enumerate() {
            let thumbnail = imageops::resize(im, thumb, thumb, imageops::FilterType::Lanczos3);
            imageops::replace(
                &mut sheet,
                &thumbnail,
                i64::from(thumb) * i64::try_from(column)?,
                i64::from(thumb) * i64::try_from(row)?,
            );
        }
        results.push(RenderComparison {
            view,
            mean_abs: metrics[0],
            changed_fraction: metrics[1],
            silhouette_xor: metrics[2],
            reference: ref_path,
            candidate: cand_path,
            diff: diff_path,
            pass: metrics
                .into_iter()
                .zip(limits)
                .all(|(value, limit)| value <= limit),
        });
    }
    sheet.save(directory.join("contact-sheet.png"))?;
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn identical_and_missing_triangle() -> Result<()> {
        let triangle = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let frame = Frame::new(&[triangle], View::Pz)?;
        let image = render(&[triangle], frame)?;
        let (_, same) = compare_pixels(&image, &image)?;
        assert!(same.into_iter().all(|value| value == 0.0));
        let empty = render(&[], frame)?;
        let (_, missing) = compare_pixels(&image, &empty)?;
        assert!(missing.into_iter().all(|value| value > 0.1));
        Ok(())
    }
    #[test]
    fn clipped_triangle_and_degenerate_edges() -> Result<()> {
        let mut image = GrayImage::from_pixel(4, 4, Luma([255]));
        fill_triangle(
            &mut image,
            &Projected {
                points: [[-10, 1], [10, 1], [0, 1]],
                depth: 0.0,
                value: 7,
            },
        )?;
        assert_eq!(
            image.as_raw(),
            &[
                255, 255, 255, 255, 7, 7, 7, 7, 255, 255, 255, 255, 255, 255, 255, 255
            ]
        );
        Ok(())
    }
}
