//! Spatial boundary integration on a sphere, without imported UV pcurves.
//!
//! For unit vectors q on the sphere, the vector area is
//! (1/2) ∮ q × dq. Solid angle is ∮ p·(q × dq)/(1+p·q),
//! modulo 4π, for any unit gauge p whose antipode misses the boundary.
//! A valid oriented spherical face selects the signed branch in (0,4π)
//! or (-4π,0). Multiple wires contribute with their own winding.
//! This diagnostic verifies curve incidence over entire rational spans,
//! uses two independently selected gauges and reports quadrature drift.

use super::super::accurate_sum;
use super::{
    spatial_boundary::{self as boundary, ClosedWire, Span, cartesian, point_derivative},
    spatial_quadrature as quadrature,
};
use crate::vector::{Point, cross, dot, norm, scale, sub};
use anyhow::{Context, Result, bail, ensure};
use occt_wasm::{OcctKernel, ShapeHandle, SurfaceData};
use serde::Serialize;

const RESIDUAL_LIMIT: f64 = 1e-11;
const FULL_SPHERE: f64 = 2.0 * std::f64::consts::TAU;
const TOLERANCES: [f64; 2] = [1e-11, 1e-13];

#[derive(Debug, Serialize)]
struct Gauge {
    direction: Point,
    denominator_margin: f64,
}

fn radial_bound(span: &[super::rational_bezier::Homogeneous]) -> Result<f64> {
    let residual = super::rational_bezier::normalized_squared_radius_bound(span)?;
    ensure!(
        residual.is_finite() && residual <= RESIDUAL_LIMIT,
        "edge not on sphere: normalized squared-radius residual {residual:e}"
    );
    Ok(residual)
}

fn gauges(spans: &[Span], radial: f64) -> Result<[Gauge; 2]> {
    let mut candidates = Vec::new();
    for x in [-1.0, 0.0, 1.0] {
        for y in [-1.0, 0.0, 1.0] {
            for z in [-1.0, 0.0, 1.0] {
                if let Some(direction) = crate::vector::unit([x, y, z], RESIDUAL_LIMIT) {
                    candidates.push(direction);
                }
            }
        }
    }
    if let Some(direction) = crate::vector::unit(
        std::array::from_fn(|i| spans.iter().map(|s| cartesian(s[0])[i]).sum()),
        RESIDUAL_LIMIT,
    ) {
        candidates.push(direction);
    }
    let norm_lower = (1.0 - radial).sqrt();
    let mut scores = candidates
        .into_iter()
        .map(|direction| {
            // Positive rational weights imply the Cartesian control hull
            // contains the entire span. This proves a margin from the gauge
            // singularity before sampling the integral.
            let min_dot = spans
                .iter()
                .flat_map(|s| s.iter())
                .map(|&h| dot(direction, cartesian(h)))
                .fold(f64::INFINITY, f64::min);
            Gauge {
                direction,
                denominator_margin: norm_lower + min_dot,
            }
        })
        .filter(|g| g.denominator_margin > 0.05)
        .collect::<Vec<_>>();
    scores.sort_by(|a, b| b.denominator_margin.total_cmp(&a.denominator_margin));
    let mut scores = scores.into_iter();
    let first = scores
        .next()
        .context("no certified nonsingular spatial gauge")?;
    let second = scores
        .find(|g| dot(g.direction, first.direction) < 0.95)
        .context("no second independent nonsingular spatial gauge")?;
    Ok([first, second])
}

fn integrand(span: &Span, t: f64, gauges: &[Gauge; 2]) -> Result<[f64; 5]> {
    let (q, dq) = point_derivative(span, t);
    let length = norm(q);
    ensure!(
        length.is_finite() && length > 0.0,
        "invalid radial projection"
    );
    let raw_cross = cross(q, dq);
    // Radial projection q/|q| cancels the derivative of the norm in
    // both cross products. No approximate UV map or tangent speed is used.
    let vector = scale(raw_cross, (length * length).recip());
    let area = gauges
        .each_ref()
        .map(|g| dot(g.direction, raw_cross) / (length * (length + dot(g.direction, q))));
    Ok([area[0], area[1], vector[0], vector[1], vector[2]])
}

#[derive(Debug, Serialize)]
struct Sample {
    signed_area: f64,
    vector_area: Point,
    volume: [f64; 4],
    unit_gauge_disagreement: f64,
    estimated_unit_quadrature_error: f64,
    estimated_volume_quadrature_error: [f64; 4],
}

fn wrap_area(area: f64, sense: f64, ambiguity: f64) -> Result<f64> {
    let positive = (sense * area).rem_euclid(FULL_SPHERE);
    ensure!(
        positive > ambiguity && positive < FULL_SPHERE - ambiguity,
        "ambiguous empty/full spherical region"
    );
    Ok(sense * positive)
}

fn sample(
    spans: &[Span],
    gauges: &[Gauge; 2],
    radius: f64,
    sense: f64,
    tolerance: f64,
    closure: f64,
) -> Result<Sample> {
    let integrals = spans
        .iter()
        .map(|span| quadrature::integrate(&|t| integrand(span, t, gauges), tolerance))
        .collect::<Result<Vec<_>>>()?;
    let sums = (0..5)
        .map(|i| accurate_sum(integrals.iter().map(|r| r.value[i])))
        .collect::<Result<Vec<_>>>()?;
    let estimated_error = accurate_sum(integrals.iter().map(|r| r.estimated_error))?;
    let min_gauge_margin = gauges
        .iter()
        .map(|g| g.denominator_margin)
        .fold(f64::INFINITY, f64::min);
    // Endpoint mismatch is measured in normalized coordinates. Treat its
    // amplification near the gauge singularity explicitly; do not silently
    // round an unresolved tiny cap into its full-sphere complement.
    let ambiguity =
        (8.0 * estimated_error + 8.0 * closure / min_gauge_margin).max(64.0 * f64::EPSILON);
    let area = wrap_area(sums[0], sense, ambiguity)?;
    let alternative = wrap_area(sums[1], sense, ambiguity)?;
    let disagreement = (area - alternative).abs();
    ensure!(
        disagreement <= ambiguity,
        "spatial gauge disagreement {disagreement:e} exceeds {ambiguity:e}"
    );
    let radius2 = radius * radius;
    Ok(Sample {
        signed_area: area * radius2,
        vector_area: [sums[2], sums[3], sums[4]].map(|x| 0.5 * radius2 * x),
        volume: [0.0; 4],
        unit_gauge_disagreement: disagreement,
        estimated_unit_quadrature_error: estimated_error,
        estimated_volume_quadrature_error: [0.0; 4],
    })
}

#[derive(Debug, Serialize)]
pub(super) struct SpatialFlux {
    center: Point,
    radius: f64,
    wires: usize,
    edges: usize,
    spans: usize,
    normalized_squared_radius_residual: f64,
    normalized_endpoint_closure_error: f64,
    gauges: [Gauge; 2],
    coarse: Sample,
    fine: Sample,
}

impl SpatialFlux {
    pub(super) fn set_references(&mut self, references: [Point; 4]) {
        for sample in [&mut self.coarse, &mut self.fine] {
            sample.volume = references.map(|reference| {
                (self.radius * sample.signed_area
                    + dot(sub(self.center, reference), sample.vector_area))
                    / 3.0
            });
            // Propagate the maximum component-wise quadrature estimate through
            // the signed flux functional. This is not an interval enclosure or
            // a bound on the geometric radial/closure residuals.
            sample.estimated_volume_quadrature_error = references.map(|reference| {
                let offset_l1: f64 = sub(self.center, reference).iter().map(|x| x.abs()).sum();
                self.radius
                    * self.radius
                    * (self.radius + 0.5 * offset_l1)
                    * sample.estimated_unit_quadrature_error
                    / 3.0
            });
        }
    }

    pub(super) const fn volumes(&self, index: usize) -> [f64; 2] {
        [self.coarse.volume[index], self.fine.volume[index]]
    }

    pub(super) const fn quadrature_error(&self, index: usize) -> [f64; 2] {
        [
            self.coarse.estimated_volume_quadrature_error[index],
            self.fine.estimated_volume_quadrature_error[index],
        ]
    }
}

fn integrate(
    center: Point,
    radius: f64,
    sense: f64,
    wires: Vec<ClosedWire>,
) -> Result<SpatialFlux> {
    ensure!(
        radius.is_finite() && radius > 0.0 && center.iter().all(|x| x.is_finite()),
        "invalid spherical carrier"
    );
    ensure!(!wires.is_empty(), "missing spherical boundary wires");
    let wire_count = wires.len();
    let edge_count = wires.iter().map(|wire| wire.edges).sum();
    let closure = accurate_sum(wires.iter().map(|wire| wire.closure_error))?;
    let spans = wires
        .into_iter()
        .flat_map(|wire| wire.spans)
        .collect::<Vec<_>>();
    let radial = spans
        .iter()
        .map(|s| radial_bound(s))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .fold(0.0, f64::max);
    let gauges = gauges(&spans, radial)?;
    let coarse = sample(&spans, &gauges, radius, sense, TOLERANCES[0], closure)?;
    let fine = sample(&spans, &gauges, radius, sense, TOLERANCES[1], closure)?;
    Ok(SpatialFlux {
        center,
        radius,
        wires: wire_count,
        edges: edge_count,
        spans: spans.len(),
        normalized_squared_radius_residual: radial,
        normalized_endpoint_closure_error: closure,
        gauges,
        coarse,
        fine,
    })
}

pub(super) fn evaluate(
    kernel: &mut OcctKernel,
    face: ShapeHandle,
    orientation: &str,
) -> Result<SpatialFlux> {
    let SurfaceData::Sphere { frame, radius } = kernel.surface_data(face)? else {
        bail!("not a spherical face");
    };
    ensure!(radius.is_finite() && radius > 0.0, "invalid sphere radius");
    let center = [frame.location.x, frame.location.y, frame.location.z];
    let sense = if (orientation == "forward") == frame.is_direct {
        1.0
    } else {
        -1.0
    };
    let wires = boundary::native_wires(kernel, face, center, radius, RESIDUAL_LIMIT)?;
    integrate(center, radius, sense, wires)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn integrate_test(
        center: Point,
        radius: f64,
        sense: f64,
        wires: Vec<Vec<Vec<Span>>>,
    ) -> Result<SpatialFlux> {
        let wires = wires
            .into_iter()
            .map(|edges| boundary::close_wire(edges, RESIDUAL_LIMIT, |_, _| Ok(true)))
            .collect::<Result<Vec<_>>>()?;
        integrate(center, radius, sense, wires)
    }

    fn latitude(z: f64, reverse: bool) -> Result<Vec<Span>> {
        let mut spans = boundary::circle_spans(
            [0.0, 0.0, z],
            (1.0 - z * z).sqrt(),
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            0.0,
            std::f64::consts::TAU,
        )?;
        if reverse {
            spans.reverse();
            for s in &mut spans {
                s.reverse();
            }
        }
        Ok(spans)
    }

    #[test]
    fn circular_caps_holes_and_complements() -> Result<()> {
        use std::f64::consts::{PI, TAU};
        let radius = 2.0;
        let cases = [
            (vec![vec![latitude(0.5, false)?]], 1.0, PI, 0.75 * PI),
            (vec![vec![latitude(0.5, true)?]], 1.0, 3.0 * PI, -0.75 * PI),
            (vec![vec![latitude(0.5, true)?]], -1.0, -PI, -0.75 * PI),
            (
                vec![vec![latitude(-0.25, false)?], vec![latitude(0.5, true)?]],
                1.0,
                TAU * 0.75,
                PI * (0.5_f64.powi(2) - 0.25_f64.powi(2)),
            ),
        ];
        for (wires, sense, area, vector_z) in cases {
            let mut flux = integrate_test([1.0, 2.0, 3.0], radius, sense, wires)?;
            flux.set_references([[1.0, 2.0, 3.0]; 4]);
            for sample in [&flux.coarse, &flux.fine] {
                assert!(
                    (sample.signed_area - radius * radius * area).abs() < 1e-10,
                    "{sample:?}"
                );
                assert!(
                    (sample.vector_area[2] - radius * radius * vector_z).abs() < 1e-10,
                    "{sample:?}"
                );
                assert!(sample.vector_area[0].abs() < 1e-10 && sample.vector_area[1].abs() < 1e-10);
                assert!((sample.volume[0] - radius.powi(3) * area / 3.0).abs() < 1e-10);
            }
        }
        Ok(())
    }

    #[test]
    fn actual_occt_small_circle_cap_and_annulus() -> Result<()> {
        use std::f64::consts::{PI, TAU};
        let mut kernel = crate::test_kernel()?;
        let sphere = kernel.make_sphere(1.0)?;
        // Rotate the poles to ±X and the seam to −Z so both trimmed
        // regions have regular closed circular wires, including the hole.
        let sphere = kernel.transform(
            sphere,
            &[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        )?;
        for (lo, hi, expected_wires) in [(0.5, 1.0, 1), (0.25, 0.5, 2)] {
            let tool = kernel.make_box_from_corners(
                -2.0,
                -2.0,
                lo,
                2.0,
                2.0,
                if hi == 1.0 { 2.0 } else { hi },
            )?;
            let region = kernel.common(sphere, tool)?;
            ensure!(kernel.is_valid(region)?, "invalid OCCT spherical region");
            let mut spherical = Vec::new();
            for face in kernel.get_sub_shape_occurrences(region, "face")? {
                if kernel.surface_type(face)? == "sphere" {
                    let orientation = kernel.shape_orientation(face)?;
                    spherical.push(evaluate(&mut kernel, face, &orientation)?);
                }
            }
            let [flux] = spherical.as_mut_slice() else {
                bail!("expected one spherical patch, got {}", spherical.len());
            };
            flux.set_references([[0.0; 3]; 4]);
            assert_eq!(flux.wires, expected_wires);
            for sample in [&flux.coarse, &flux.fine] {
                assert!(
                    (sample.signed_area - TAU * (hi - lo)).abs() < 1e-11,
                    "{sample:?}"
                );
                assert!(
                    (sample.vector_area[2] - PI * (hi * hi - lo * lo)).abs() < 1e-11,
                    "{sample:?}"
                );
                assert!(sample.vector_area[0].abs() < 1e-11 && sample.vector_area[1].abs() < 1e-11);
            }
        }
        Ok(())
    }

    #[test]
    fn rejects_an_off_sphere_rational_boundary() -> Result<()> {
        let mut edge = latitude(0.25, false)?;
        edge[0][1][2] += 1e-5;
        assert!(integrate_test([0.0; 3], 1.0, 1.0, vec![vec![edge]]).is_err());
        Ok(())
    }
}
