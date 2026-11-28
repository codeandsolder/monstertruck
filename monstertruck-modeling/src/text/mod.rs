//! Font outline ingestion and conversion to `monstertruck` curves and wires.
//!
//! Requires the `font` feature flag. Parses TrueType/OpenType glyph outlines
//! via [`skrifa`] and maps line, quadratic Bezier, and cubic Bezier
//! segments into [`Curve`](crate::Curve) edges assembled into closed
//! [`Wire`](crate::Wire)s suitable for planar face construction.

use crate::{
    Curve, Cut, MetricSpace, ParametricCurve, Point3, Result, Vertex, builder, errors::Error,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use monstertruck_geometry::prelude::{BsplineCurve, KnotVector};
use skrifa::{DrawMemory, PathSink, ScalerFont};
pub use skrifa::{Font, GlyphId};

type Edge = monstertruck_topology::Edge<Point3, Curve>;
type Wire = monstertruck_topology::Wire<Point3, Curve>;

/// Options for glyph/text profile extraction.
#[derive(Debug, Clone)]
pub struct TextOptions {
    /// Uniform scale applied to font units. Default: `1.0 / units_per_em`.
    pub scale: Option<f64>,
    /// Whether to flip the Y axis (font coordinates have Y-up, but many CAD
    /// systems use Y-down or want the baseline at `y = 0` growing upward).
    /// Default: `true` (flip font Y-up coordinates to Y-down).
    pub y_flip: bool,
    /// Z coordinate for the planar profile. Default: `0.0`.
    pub z: f64,
    /// Tolerance for closing contours: if the distance between the last point
    /// and the contour start is below this, they are snapped together.
    /// Default: `1e-7`.
    pub closure_tolerance: f64,
}

impl Default for TextOptions {
    fn default() -> Self {
        Self {
            scale: None,
            y_flip: true,
            z: 0.0,
            closure_tolerance: DEFAULT_NUMERICAL_TOLERANCES.text_closure,
        }
    }
}

/// A single outline segment extracted from a font glyph.
#[derive(Debug, Clone)]
enum Segment {
    /// Straight line to `(x, y)`.
    Line(f64, f64),
    /// Quadratic Bezier with control point `(cx, cy)` to endpoint `(x, y)`.
    Quad(f64, f64, f64, f64),
    /// Cubic Bezier with control points `(c1x, c1y)`, `(c2x, c2y)` to `(x, y)`.
    Cubic(f64, f64, f64, f64, f64, f64),
}

/// Collects glyph outline contours from [`PathSink`] callbacks.
struct ContourCollector {
    contours: Vec<(f64, f64, Vec<Segment>)>,
    current_start: Option<(f64, f64)>,
    current_segments: Vec<Segment>,
}

impl ContourCollector {
    const fn new() -> Self {
        Self {
            contours: Vec::new(),
            current_start: None,
            current_segments: Vec::new(),
        }
    }
}

impl PathSink for ContourCollector {
    fn move_to(&mut self, x: f32, y: f32) {
        self.close();
        self.current_start = Some((f64::from(x), f64::from(y)));
    }

    fn line_to(&mut self, x: f32, y: f32) {
        self.current_segments
            .push(Segment::Line(f64::from(x), f64::from(y)));
    }

    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        self.current_segments.push(Segment::Quad(
            f64::from(x1),
            f64::from(y1),
            f64::from(x),
            f64::from(y),
        ));
    }

    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        self.current_segments.push(Segment::Cubic(
            f64::from(x1),
            f64::from(y1),
            f64::from(x2),
            f64::from(y2),
            f64::from(x),
            f64::from(y),
        ));
    }

    fn close(&mut self) {
        if let Some((sx, sy)) = self.current_start.take() {
            let segs = std::mem::take(&mut self.current_segments);
            if !segs.is_empty() {
                self.contours.push((sx, sy, segs));
            }
        }
    }
}

/// Placement in the output coordinate system, evaluated uniformly for
/// every endpoint and control point of a contour.
#[derive(Clone, Copy)]
struct Placement {
    scale: f64,
    y_flip: bool,
    z: f64,
    closure_tolerance: f64,
    offset: [f64; 2],
}

impl Placement {
    fn new(font: &Font, opts: &TextOptions) -> Result<Self> {
        let scale = opts
            .scale
            .unwrap_or_else(|| 1.0 / f64::from(font.units_per_em()));
        if !scale.is_finite()
            || scale <= 0.0
            || !opts.z.is_finite()
            || !opts.closure_tolerance.is_finite()
            || opts.closure_tolerance < 0.0
        {
            return Err(Error::InvalidTextOptions);
        }
        Ok(Self {
            scale,
            y_flip: opts.y_flip,
            z: opts.z,
            closure_tolerance: opts.closure_tolerance,
            offset: [0.0; 2],
        })
    }

    fn point(self, x: f64, y: f64) -> Point3 {
        let y = self.scale.mul_add(y, self.offset[1]);
        Point3::new(
            self.scale.mul_add(x, self.offset[0]),
            if self.y_flip { -y } else { y },
            self.z,
        )
    }

    fn positioned(self, position: [f64; 2]) -> Result<Self> {
        let offset = position.map(|value| value * self.scale);
        if !offset.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidTextOptions);
        }
        Ok(Self { offset, ..self })
    }
}

struct PlacedSegment {
    end: Point3,
    controls: Vec<Point3>,
}

impl Segment {
    fn place(&self, placement: Placement) -> PlacedSegment {
        let (controls, (x, y)) = match *self {
            Self::Line(x, y) => (vec![], (x, y)),
            Self::Quad(cx, cy, x, y) => (vec![placement.point(cx, cy)], (x, y)),
            Self::Cubic(cx, cy, dx, dy, x, y) => (
                vec![placement.point(cx, cy), placement.point(dx, dy)],
                (x, y),
            ),
        };
        PlacedSegment {
            end: placement.point(x, y),
            controls,
        }
    }
}

fn finite_point(point: Point3) -> bool {
    [point.x, point.y, point.z].into_iter().all(f64::is_finite)
}

/// A closed, nonconstant Bezier segment needs two distinct topological edges.
/// Subdivide the polynomial exactly rather than introducing a line or an
/// invalid same-vertex edge.
fn append_closed_bezier(
    edges: &mut Vec<Edge>,
    vertex: &Vertex,
    mut controls: Vec<Point3>,
) -> Result<()> {
    let point = vertex.point();
    controls.insert(0, point);
    controls.push(point);
    let mut curve =
        BsplineCurve::new_unchecked(KnotVector::bezier_knot(controls.len() - 1), controls);
    // A nonconstant polynomial of degree at most three cannot equal its
    // starting point at all three interior probes and both endpoints.
    let (parameter, midpoint) = [0.5, 0.25, 0.75]
        .into_iter()
        .map(|t| (t, curve.evaluate(t)))
        .find(|(_, p)| *p != point && finite_point(*p))
        .ok_or_else(empty_outline)?;
    let middle = builder::vertex(midpoint);
    let right = curve.cut(parameter);
    edges.push(Edge::new(vertex, &middle, Curve::BsplineCurve(curve))?);
    edges.push(Edge::new(&middle, vertex, Curve::BsplineCurve(right))?);
    Ok(())
}

/// Converts a collected contour into a closed wire.
fn contour_to_wire(
    start_x: f64,
    start_y: f64,
    segments: &[Segment],
    placement: Placement,
) -> Result<Wire> {
    let first_pt = placement.point(start_x, start_y);
    if !finite_point(first_pt) {
        return Err(Error::InvalidTextOptions);
    }
    let first_vertex = builder::vertex(first_pt);
    let mut edges = Vec::with_capacity(segments.len());
    let mut prev_vertex = first_vertex.clone();
    let mut previous = first_pt;
    let mut placed = segments
        .iter()
        .map(|seg| seg.place(placement))
        .filter(|seg| {
            let constant = seg.end == previous && seg.controls.iter().all(|cp| *cp == previous);
            previous = seg.end;
            !constant
        })
        .peekable();
    while let Some(segment) = placed.next() {
        let PlacedSegment { end, controls } = segment;
        if !finite_point(end) || !controls.iter().copied().all(finite_point) {
            return Err(Error::InvalidTextOptions);
        }
        let next =
            if placed.peek().is_none() && first_pt.distance(end) <= placement.closure_tolerance {
                first_vertex.clone()
            } else {
                builder::vertex(end)
            };
        if prev_vertex == next {
            if controls.iter().any(|cp| *cp != next.point()) {
                append_closed_bezier(&mut edges, &next, controls)?;
            }
        } else {
            let edge: Edge = if controls.is_empty() {
                builder::line(&prev_vertex, &next)
            } else {
                builder::bezier(&prev_vertex, &next, controls)
            };
            edges.push(edge);
        }
        prev_vertex = next;
    }
    if prev_vertex.id() != first_vertex.id() {
        edges.push(builder::line(&prev_vertex, &first_vertex));
    }
    if edges.is_empty() {
        return Err(empty_outline());
    }
    Ok(edges.into())
}

const fn empty_outline() -> Error {
    Error::FromTopology(monstertruck_topology::errors::Error::EmptyWire)
}

/// Skips degenerate contours where all segments collapse to a single point.
fn is_degenerate_contour(start_x: f64, start_y: f64, segments: &[Segment]) -> bool {
    let same = |x: f64, y: f64| {
        (x - start_x).abs() < DEFAULT_NUMERICAL_TOLERANCES.structural
            && (y - start_y).abs() < DEFAULT_NUMERICAL_TOLERANCES.structural
    };
    segments.iter().all(|seg| match *seg {
        Segment::Line(x, y) => same(x, y),
        Segment::Quad(cx, cy, x, y) => same(cx, cy) && same(x, y),
        Segment::Cubic(cx, cy, dx, dy, x, y) => same(cx, cy) && same(dx, dy) && same(x, y),
    })
}

/// A glyph and its baseline-relative origin in unscaled font units.
/// Callers may supply kerning, shaping, ligatures and vertical offsets.
#[derive(Debug, Clone, Copy)]
pub struct PositionedGlyph {
    /// Glyph ID in the selected font instance.
    pub glyph_id: GlyphId,
    /// X and Y offset in the font's coordinate system (Y-up).
    pub position: [f64; 2],
}

struct PlacedContour {
    start_x: f64,
    start_y: f64,
    segments: Vec<Segment>,
    placement: Placement,
}

impl PlacedContour {
    fn to_wire(&self) -> Result<Wire> {
        contour_to_wire(self.start_x, self.start_y, &self.segments, self.placement)
    }
}

/// Converts a single glyph outline to closed wires, including interior holes.
/// Uses the selected font instance and unhinted font-unit geometry.
///
/// # Errors
/// Returns an error for an invalid placement, missing/degenerate outline,
/// or a malformed glyph.
pub fn glyph_profile(font: &Font, glyph_id: GlyphId, opts: &TextOptions) -> Result<Vec<Wire>> {
    let wires = positioned_glyph_profiles(
        font,
        [PositionedGlyph {
            glyph_id,
            position: [0.0; 2],
        }],
        opts,
    )?;
    if wires.is_empty() {
        Err(empty_outline())
    } else {
        Ok(wires)
    }
}

/// Converts an externally positioned glyph run to planar wires.
///
/// The font's selected collection face and variation instance are respected.
/// Every coordinate is transformed by the same placement. Glyphs without
/// contours (such as spaces) contribute no wires.
///
/// # Errors
/// Returns an error for invalid placement or a malformed outline.
pub fn positioned_glyph_profiles(
    font: &Font,
    glyphs: impl IntoIterator<Item = PositionedGlyph>,
    opts: &TextOptions,
) -> Result<Vec<Wire>> {
    let placement = Placement::new(font, opts)?;
    let scaler = ScalerFont::unhinted(font.clone(), None);
    let mut contours = Vec::new();
    let mut memory = Vec::new();
    for glyph in glyphs {
        let placement = placement.positioned(glyph.position)?;
        if glyph.glyph_id.to_u32() >= font.num_glyphs() {
            return Err(Error::InvalidGlyphOutline(glyph.glyph_id.to_u32()));
        }
        let Some(outline) = scaler.outline(glyph.glyph_id) else {
            continue;
        };
        let mut collector = ContourCollector::new();
        memory.resize(outline.required_memory_size(), 0);
        outline
            .draw(DrawMemory::Buffer(&mut memory), &mut collector)
            .map_err(|_| Error::InvalidGlyphOutline(glyph.glyph_id.to_u32()))?;
        collector.close();
        contours.extend(
            collector
                .contours
                .into_iter()
                .filter(|(sx, sy, segs)| !is_degenerate_contour(*sx, *sy, segs))
                .map(|(start_x, start_y, segments)| PlacedContour {
                    start_x,
                    start_y,
                    segments,
                    placement,
                }),
        );
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        use rayon::prelude::*;
        contours.par_iter().map(PlacedContour::to_wire).collect()
    }
    #[cfg(target_arch = "wasm32")]
    {
        contours.iter().map(PlacedContour::to_wire).collect()
    }
}

/// Converts Unicode characters to glyph outlines using their horizontal advances.
///
/// This is a simple cmap/advance layout, not a shaping or kerning engine.
/// Use [`positioned_glyph_profiles`] for externally shaped text.
/// Unmapped characters are skipped; spaces retain their advance.
///
/// # Errors
/// Returns an error for invalid placement or a malformed outline.
pub fn text_profile(font: &Font, text: &str, opts: &TextOptions) -> Result<Vec<Wire>> {
    let charmap = font.charmap();
    let metrics = font.glyph_metrics();
    let mut cursor = 0.0;
    let glyphs = text
        .chars()
        .filter_map(|ch| charmap.map_unicode(ch))
        .map(|glyph_id| {
            let glyph = PositionedGlyph {
                glyph_id,
                position: [cursor, 0.0],
            };
            cursor += f64::from(metrics.h_advance(glyph_id));
            glyph
        });
    positioned_glyph_profiles(font, glyphs, opts)
}

#[cfg(test)]
mod tests;
