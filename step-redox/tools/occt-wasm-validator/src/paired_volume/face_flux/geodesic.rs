//! Analytic flux for a convex minor-geodesic spherical polygon. This is a
//! diagnostic alternative to the UV quadrature, not an admission criterion.
//! Every rational edge is checked over the entire finite trim using
//! homogeneous Bernstein bounds (computed in floating point).

use super::{
    array,
    rational_bezier::{self, Homogeneous},
};
use anyhow::{Context, Result, bail, ensure};
use occt_wasm::{CurveGeometry, OcctKernel, ShapeHandle, SurfaceData, Vec3};
use serde::Serialize;

use crate::vector::{Point, cross, dot, norm, scale, sub};

const RESIDUAL_LIMIT: f64 = 1.0e-11;

fn unit(a: Point) -> Result<Point> {
    crate::vector::unit(a, RESIDUAL_LIMIT).context("degenerate direction")
}
fn xyz(v: Vec3) -> Point { [v.x, v.y, v.z] }
fn spatial(h: Homogeneous) -> Point { [h[0], h[1], h[2]] }

#[derive(Debug, Serialize)]
struct Arc {
    start: Point,
    end: Point,
    normal: Point,
    angle: f64,
    plane_residual: f64,
    squared_radius_residual: f64,
}

fn check_spline(
    spans: &[Vec<Homogeneous>],
    start: Point,
    end: Point,
    normal: Point,
) -> Result<(f64, f64)> {
    let mut plane = 0.0_f64;
    let mut radial = 0.0_f64;
    for poles in spans {
        let min_weight = poles.iter().map(|p| p[3]).fold(f64::INFINITY, f64::min);
        ensure!(min_weight > 0.0, "nonpositive refined weight");
        for &pole in poles {
            let point = scale(spatial(pole), 1.0 / pole[3]);
            plane = plane.max(dot(point, normal).abs());
            ensure!(
                dot(cross(start, point), normal) >= -RESIDUAL_LIMIT
                    && dot(cross(point, end), normal) >= -RESIDUAL_LIMIT,
                "spline control outside minor arc wedge"
            );
        }
        let coefficients = rational_bezier::product(poles, poles, |a, b| {
            dot(spatial(a), spatial(b)) - a[3] * b[3]
        })?;
        radial = radial
            .max(coefficients.into_iter().map(f64::abs).fold(0.0, f64::max) / min_weight.powi(2));
        let differences = poles
            .windows(2)
            .map(|p| std::array::from_fn(|i| p[1][i] - p[0][i]))
            .collect::<Vec<_>>();
        let angular = rational_bezier::product(poles, &differences, |a, b| {
            dot(normal, cross(spatial(a), spatial(b)))
        })?;
        ensure!(
            angular.iter().all(|&x| x > 0.0),
            "angular progression is not strictly monotone"
        );
    }
    ensure!(
        plane <= RESIDUAL_LIMIT && radial <= RESIDUAL_LIMIT,
        "not a great-circle spline: plane={plane:e}, squared_radius={radial:e}"
    );
    Ok((plane, radial))
}

fn arc(kernel: &mut OcctKernel, edge: ShapeHandle, center: Point, radius: f64) -> Result<Arc> {
    let data = kernel.curve_data(edge)?;
    let first = data.first;
    let last = data.last;
    ensure!(
        first.is_finite() && last.is_finite() && first < last,
        "invalid finite trim"
    );
    let start = scale(
        sub(
            array(kernel.curve_point_at_param(edge, first)?, "curve endpoint")?,
            center,
        ),
        1.0 / radius,
    );
    let end = scale(
        sub(
            array(kernel.curve_point_at_param(edge, last)?, "curve endpoint")?,
            center,
        ),
        1.0 / radius,
    );
    ensure!(
        (norm(start) - 1.0).abs() <= RESIDUAL_LIMIT && (norm(end) - 1.0).abs() <= RESIDUAL_LIMIT,
        "endpoint not on sphere"
    );
    let start = unit(start)?;
    let end = unit(end)?;
    let cr = cross(start, end);
    let normal = unit(cr)?;
    let angle = norm(cr).atan2(dot(start, end));
    let (plane_residual, squared_radius_residual) = match data.geometry {
        CurveGeometry::Circle {
            center: c,
            radius: r,
            axis,
            ..
        } => {
            let offset = norm(sub(xyz(c), center)) / radius;
            let dr = (r / radius - 1.0).abs();
            let normal_error = norm(cross(unit(xyz(axis))?, normal));
            ensure!(
                offset <= RESIDUAL_LIMIT
                    && dr <= RESIDUAL_LIMIT
                    && normal_error <= RESIDUAL_LIMIT
                    && (last - first - angle).abs() <= RESIDUAL_LIMIT,
                "circle does not trace the minor great-circle arc"
            );
            (offset + normal_error, (offset + dr) * (2.0 + offset + dr))
        }
        CurveGeometry::BSpline(ref spline) => {
            let spans = rational_bezier::spans(spline, center, radius, first, last)?;
            check_spline(&spans, start, end, normal)?
        }
        _ => bail!("unsupported geodesic edge carrier"),
    };
    let mut arc = Arc {
        start,
        end,
        normal,
        angle,
        plane_residual,
        squared_radius_residual,
    };
    match kernel.shape_orientation(edge)?.as_str() {
        "forward" => (),
        "reversed" => {
            std::mem::swap(&mut arc.start, &mut arc.end);
            arc.normal = scale(arc.normal, -1.0);
        }
        _ => bail!("non-oriented geodesic edge"),
    }
    Ok(arc)
}

#[derive(Debug, Serialize)]
pub(super) struct GeodesicFlux {
    center: Point,
    radius: f64,
    arcs: Vec<Arc>,
    signed_area: f64,
    vector_area: Point,
    volume: [f64; 4],
}

impl GeodesicFlux {
    pub(super) fn set_references(&mut self, references: [Point; 4]) {
        self.volume = references.map(|r| {
            (self.radius * self.signed_area + dot(sub(self.center, r), self.vector_area)) / 3.0
        });
    }
}

fn polygon(center: Point, radius: f64, mut arcs: Vec<Arc>, sense: f64) -> Result<GeodesicFlux> {
    ensure!(arcs.len() >= 3, "too few geodesic edges");
    let mut ordered = vec![arcs.pop().context("missing edge")?];
    while !arcs.is_empty() {
        let end = ordered.last().context("missing previous edge")?.end;
        let mut candidates = arcs
            .iter()
            .enumerate()
            .filter(|(_, a)| norm(sub(a.start, end)) <= RESIDUAL_LIMIT);
        let (i, _) = candidates.next().context("unclosed directed boundary")?;
        ensure!(candidates.next().is_none(), "ambiguous directed boundary");
        ordered.push(arcs.swap_remove(i));
    }
    ensure!(
        norm(sub(
            ordered[0].start,
            ordered.last().context("missing last edge")?.end
        )) <= RESIDUAL_LIMIT,
        "unclosed directed boundary"
    );
    let vertices = ordered.iter().map(|a| a.start).collect::<Vec<_>>();
    let hemisphere = unit(std::array::from_fn(|i| vertices.iter().map(|p| p[i]).sum()))?;
    ensure!(
        vertices
            .iter()
            .all(|&p| dot(p, hemisphere) > RESIDUAL_LIMIT),
        "not contained in a minor hemisphere"
    );
    // The oriented half-spaces must all enclose the same convex region.
    // A reversed face must have reversed boundary winding, not the complement.
    for arc in &ordered {
        ensure!(
            vertices
                .iter()
                .all(|&p| sense * dot(arc.normal, p) >= -RESIDUAL_LIMIT),
            "nonconvex or complementary spherical boundary"
        );
    }
    let solid_angle: f64 = ordered
        .iter()
        .map(|a| {
            2.0 * dot(hemisphere, cross(a.start, a.end)).atan2(
                1.0 + dot(hemisphere, a.start) + dot(a.start, a.end) + dot(a.end, hemisphere),
            )
        })
        .sum();
    ensure!(
        sense * solid_angle > RESIDUAL_LIMIT && solid_angle.abs() < std::f64::consts::TAU,
        "invalid oriented solid angle"
    );
    let vector_area = std::array::from_fn(|i| {
        0.5 * radius * radius * ordered.iter().map(|a| a.angle * a.normal[i]).sum::<f64>()
    });
    Ok(GeodesicFlux {
        center,
        radius,
        arcs: ordered,
        signed_area: solid_angle * radius * radius,
        vector_area,
        volume: [0.0; 4],
    })
}

pub(super) fn evaluate(
    kernel: &mut OcctKernel,
    face: ShapeHandle,
    orientation: &str,
) -> Result<GeodesicFlux> {
    let SurfaceData::Sphere { frame, radius } = kernel.surface_data(face)? else {
        bail!("not a spherical face");
    };
    ensure!(radius.is_finite() && radius > 0.0, "invalid sphere radius");
    let center = xyz(frame.location);
    let wires = kernel.get_sub_shape_occurrences(face, "wire")?;
    let [wire] = wires.as_slice() else {
        bail!("not a single-wire spherical face");
    };
    let edges = kernel.get_sub_shape_occurrences(*wire, "edge")?;
    let arcs = edges
        .into_iter()
        .map(|e| arc(kernel, e, center, radius))
        .collect::<Result<Vec<_>>>()?;
    let sense = if (orientation == "forward") == frame.is_direct {
        1.0
    } else {
        -1.0
    };
    polygon(center, radius, arcs, sense)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_quarter_circle_checks_entire_span() -> Result<()> {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let start = [1.0, 0.0, 0.0];
        let end = [0.0, 1.0, 0.0];
        let normal = [0.0, 0.0, 1.0];
        let good = vec![[1.0, 0.0, 0.0, 1.0], [w, w, 0.0, w], [0.0, 1.0, 0.0, 1.0]];
        let (plane, radial) = check_spline(&[good.clone()], start, end, normal)?;
        assert!(plane < 1e-15 && radial < 1e-15);
        let mut bad = good;
        bad[1][0] += 1e-5;
        assert!(check_spline(&[bad], start, end, normal).is_err());
        Ok(())
    }

    #[test]
    fn octant_flux_and_reversed_winding() -> Result<()> {
        let points = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let arcs = |reverse: bool| {
            (0..3)
                .map(|i| {
                    let mut start = points[i];
                    let mut end = points[(i + 1) % 3];
                    if reverse {
                        std::mem::swap(&mut start, &mut end);
                    }
                    Arc {
                        start,
                        end,
                        normal: cross(start, end),
                        angle: std::f64::consts::FRAC_PI_2,
                        plane_residual: 0.0,
                        squared_radius_residual: 0.0,
                    }
                })
                .collect()
        };
        let mut flux = polygon([0.0; 3], 1.0, arcs(false), 1.0)?;
        flux.set_references([[0.0; 3]; 4]);
        assert!((flux.volume[0] - std::f64::consts::PI / 6.0).abs() < 1e-15);
        assert!(
            flux.vector_area
                .iter()
                .all(|v| (v - std::f64::consts::FRAC_PI_4).abs() < 1e-15)
        );
        assert!(polygon([0.0; 3], 1.0, arcs(false), -1.0).is_err());
        let mut reverse = polygon([0.0; 3], 1.0, arcs(true), -1.0)?;
        reverse.set_references([[0.0; 3]; 4]);
        assert!((reverse.volume[0] + flux.volume[0]).abs() < 1e-15);
        Ok(())
    }
}
