//! Independent spatial flux on elementary carriers. Uses native OCCT
//! topology and actual trimmed 3-D curves, never reconstructed CAD IR,
//! UV pcurves, meshing, or the signed-face surface quadrature.
//!
//! Small radial/planar residuals and endpoint connectors are reported,
//! not silently treated as exact geometry. The resulting paired comparison
//! is a diagnostic, not an interval proof or a substitute for the raw gate.

use super::{
    rational_bezier::{self, Homogeneous},
    spatial_boundary::{self, Span, cartesian, point_derivative},
    spatial_quadrature as quadrature, spatial_sphere,
};
use crate::{
    paired_volume::{accurate_sum, references},
    vector::{Point, add, cross, dot, norm, scale, sub, unit, xyz},
};
use anyhow::{Context, Result, bail, ensure};
use occt_wasm::{OcctKernel, ShapeHandle, SurfaceData};
use serde::Serialize;

const QUADRATURE_EPS: [f64; 2] = [1e-11, 1e-13];

#[derive(Clone, Copy)]
enum Carrier {
    Plane {
        normal: Point,
        origin: Point,
        length: f64,
    },
    Cylinder {
        axis: Point,
        origin: Point,
        radius: f64,
    },
}

impl Carrier {
    fn frame(self) -> (Point, f64) {
        match self {
            Self::Plane { origin, length, .. } => (origin, length),
            Self::Cylinder { origin, radius, .. } => (origin, radius),
        }
    }

    fn incidence(self, span: &[Homogeneous]) -> Result<f64> {
        let (_, length) = self.frame();
        match self {
            Self::Plane { normal, .. } => {
                // Positive rational basis functions give a convex hull
                // enclosure of the signed plane distances.
                Ok(length
                    * span
                        .iter()
                        .map(|&h| dot(cartesian(h), normal).abs())
                        .fold(0.0, f64::max))
            }
            Self::Cylinder { axis, .. } => {
                let radial = span
                    .iter()
                    .map(|h| {
                        let q = rational_bezier::spatial(*h);
                        let r = sub(q, scale(axis, dot(q, axis)));
                        [r[0], r[1], r[2], h[3]]
                    })
                    .collect::<Vec<_>>();
                let squared = rational_bezier::normalized_squared_radius_bound(&radial)?;
                ensure!(
                    squared < 1.0,
                    "radial projection is not separated from the axis"
                );
                Ok(length * squared / (1.0 + (1.0 - squared).sqrt()))
            }
        }
    }

    /// Return signed area and vector-area integrands, with a separate
    /// angular winding term on a cylinder. The parameterized surface is
    /// X = O + R(n(theta) + z*a), with oriented area R²*n*dtheta*dz.
    /// Green's theorem gives int f(theta)dtheta dz = -int z*f(theta)dtheta.
    /// The formula remains valid across the periodic seam and for holes.
    fn integrand(self, span: &Span, t: f64) -> Result<[f64; 5]> {
        let (q, dq) = point_derivative(span, t);
        match self {
            Self::Plane { normal, .. } => {
                let area = 0.5 * dot(normal, cross(q, dq));
                let vector = scale(normal, area);
                Ok([area, vector[0], vector[1], vector[2], 0.0])
            }
            Self::Cylinder { axis, .. } => {
                let height = dot(q, axis);
                let radial = sub(q, scale(axis, height));
                let length = norm(radial);
                ensure!(
                    length.is_finite() && length > 0.0,
                    "singular cylindrical boundary"
                );
                let angular = dot(axis, cross(radial, dq)) / (length * length);
                let area = -height * angular;
                let vector = scale(radial, area / length);
                Ok([area, vector[0], vector[1], vector[2], angular])
            }
        }
    }

    fn sample(
        self,
        spans: &[Span],
        sense: f64,
        epsilon: f64,
        closure: f64,
        references: [Point; 4],
    ) -> Result<Sample> {
        let terms = spans
            .iter()
            .map(|span| {
                quadrature::integrate(&|t| self.integrand(span, t), epsilon / spans.len() as f64)
            })
            .collect::<Result<Vec<_>>>()?;
        let sums = (0..5)
            .map(|i| accurate_sum(terms.iter().map(|r| r.value[i])))
            .collect::<Result<Vec<_>>>()?;
        let error = accurate_sum(terms.iter().map(|r| r.estimated_error))?;
        let (origin, length) = self.frame();
        ensure!(
            sense * sums[0] > 0.0,
            "spatial wire orientation disagrees with native face"
        );
        if matches!(self, Self::Cylinder { .. }) {
            ensure!(
                sums[4].abs() <= (8.0 * (error + closure)).max(128.0 * f64::EPSILON),
                "unbalanced cylindrical angular winding {:e}",
                sums[4]
            );
        }
        let signed_area = sums[0] * length * length;
        let vector_area = scale([sums[1], sums[2], sums[3]], length * length);
        let intrinsic_length = match self {
            Self::Plane { .. } => 0.0,
            Self::Cylinder { radius, .. } => radius,
        };
        let volume = references.map(|reference| {
            (intrinsic_length * signed_area + dot(sub(origin, reference), vector_area)) / 3.0
        });
        let volume_error = references.map(|reference| {
            let offset_l1 = sub(origin, reference).iter().map(|x| x.abs()).sum::<f64>();
            (intrinsic_length + offset_l1) * length * length * error / 3.0
        });
        Ok(Sample {
            signed_area,
            vector_area,
            angular_winding: sums[4],
            volume,
            volume_error,
        })
    }
}

#[derive(Debug, Serialize)]
struct Sample {
    signed_area: f64,
    vector_area: Point,
    angular_winding: f64,
    volume: [f64; 4],
    volume_error: [f64; 4],
}

#[derive(Debug, Serialize)]
struct ElementaryFlux {
    carrier_origin: Point,
    normalization_length_mm: f64,
    wires: usize,
    edges: usize,
    spans: usize,
    incidence_upper_mm: f64,
    endpoint_closure_sum_mm: f64,
    coarse: Sample,
    fine: Sample,
}

fn elementary(
    kernel: &mut OcctKernel,
    face: ShapeHandle,
    carrier: Carrier,
    sense: f64,
    boundary_tolerance: f64,
    references: [Point; 4],
) -> Result<ElementaryFlux> {
    let (origin, length) = carrier.frame();
    let wires =
        spatial_boundary::native_wires(kernel, face, origin, length, boundary_tolerance / length)?;
    let wire_count = wires.len();
    let edges = wires.iter().map(|w| w.edges).sum();
    let closure = accurate_sum(wires.iter().map(|w| w.closure_error))?;
    let spans = wires.into_iter().flat_map(|w| w.spans).collect::<Vec<_>>();
    ensure!(!spans.is_empty(), "empty spatial boundary");
    let incidence = spans
        .iter()
        .map(|span| carrier.incidence(span))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .fold(0.0, f64::max);
    ensure!(
        incidence.is_finite() && incidence <= boundary_tolerance,
        "carrier incidence upper bound {incidence:e} mm exceeds {boundary_tolerance:e} mm"
    );
    let [coarse, fine] =
        QUADRATURE_EPS.map(|epsilon| carrier.sample(&spans, sense, epsilon, closure, references));
    Ok(ElementaryFlux {
        carrier_origin: origin,
        normalization_length_mm: length,
        wires: wire_count,
        edges,
        spans: spans.len(),
        incidence_upper_mm: incidence,
        endpoint_closure_sum_mm: closure * length,
        coarse: coarse?,
        fine: fine?,
    })
}

#[derive(Debug, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
enum Flux {
    Elementary {
        flux: Box<ElementaryFlux>,
    },
    Sphere {
        flux: Box<spatial_sphere::SpatialFlux>,
    },
    Unsupported {
        reason: String,
    },
}

impl Flux {
    fn terms(&self, index: usize) -> Option<([f64; 2], [f64; 2])> {
        match self {
            Self::Elementary { flux } => Some((
                [flux.coarse.volume[index], flux.fine.volume[index]],
                [
                    flux.coarse.volume_error[index],
                    flux.fine.volume_error[index],
                ],
            )),
            Self::Sphere { flux } => Some((flux.volumes(index), flux.quadrature_error(index))),
            Self::Unsupported { .. } => None,
        }
    }
}

fn face_flux(
    kernel: &mut OcctKernel,
    face: ShapeHandle,
    orientation: &str,
    boundary_tolerance: f64,
    references: [Point; 4],
) -> Result<Flux> {
    let data = kernel.surface_data(face)?;
    let frame = match data {
        SurfaceData::Plane(frame) | SurfaceData::Cylinder { frame, .. } => frame,
        SurfaceData::Sphere { .. } => {
            let mut flux = spatial_sphere::evaluate(kernel, face, orientation)?;
            flux.set_references(references);
            return Ok(Flux::Sphere {
                flux: Box::new(flux),
            });
        }
        _ => bail!("unsupported elementary spatial carrier"),
    };
    let axis = unit(xyz(frame.axis), 1e-14).context("invalid carrier axis")?;
    let location = xyz(frame.location);
    let offset = sub(references[0], location);
    let sense = if (orientation == "forward") == frame.is_direct {
        1.0
    } else {
        -1.0
    };
    let carrier = match data {
        SurfaceData::Plane(_) => Carrier::Plane {
            normal: axis,
            // Move the origin within the plane to avoid subtracting large
            // unrelated coordinates on remote placement frames.
            origin: sub(references[0], scale(axis, dot(offset, axis))),
            length: 1.0,
        },
        SurfaceData::Cylinder { radius, .. } => {
            ensure!(
                radius.is_finite() && radius > boundary_tolerance,
                "invalid cylinder radius"
            );
            Carrier::Cylinder {
                axis,
                origin: add(location, scale(axis, dot(offset, axis))),
                radius,
            }
        }
        _ => bail!("unsupported carrier"),
    };
    Ok(Flux::Elementary {
        flux: Box::new(elementary(
            kernel,
            face,
            carrier,
            sense,
            boundary_tolerance,
            references,
        )?),
    })
}

#[derive(Debug, Serialize)]
struct Face {
    shell_occurrence: usize,
    face_occurrence: usize,
    surface_kind: String,
    orientation: String,
    #[serde(flatten)]
    flux: Flux,
}

fn faces(
    kernel: &mut OcctKernel,
    shape: ShapeHandle,
    boundary_tolerance: f64,
    references: [Point; 4],
) -> Result<Vec<Face>> {
    let solids = kernel.get_sub_shape_occurrences(shape, "solid")?;
    let [solid] = solids.as_slice() else {
        bail!("expected one native solid occurrence");
    };
    let mut faces = Vec::new();
    for (shell_occurrence, shell) in kernel
        .get_sub_shape_occurrences(*solid, "shell")?
        .into_iter()
        .enumerate()
    {
        for (face_occurrence, face) in kernel
            .get_sub_shape_occurrences(shell, "face")?
            .into_iter()
            .enumerate()
        {
            let orientation = kernel.shape_orientation(face)?;
            if matches!(orientation.as_str(), "internal" | "external") {
                continue;
            }
            ensure!(
                matches!(orientation.as_str(), "forward" | "reversed"),
                "unknown face orientation"
            );
            let surface_kind = kernel.surface_type(face)?;
            let flux = face_flux(kernel, face, &orientation, boundary_tolerance, references)
                .unwrap_or_else(|error| Flux::Unsupported {
                    reason: format!("{error:#}"),
                });
            faces.push(Face {
                shell_occurrence,
                face_occurrence,
                surface_kind,
                orientation,
                flux,
            });
        }
    }
    ensure!(!faces.is_empty(), "missing oriented faces");
    Ok(faces)
}

#[derive(Debug, Serialize)]
struct ComparisonSample {
    reference: Point,
    volumes: [f64; 2],
    delta: f64,
    coarse_delta: f64,
    quadrature_drift: [f64; 2],
    estimated_quadrature_error: f64,
}

fn compare_sample(
    a: &[Face],
    b: &[Face],
    index: usize,
    reference: Point,
) -> Result<ComparisonSample> {
    let terms = |faces: &[Face]| {
        faces
            .iter()
            .map(|face| face.flux.terms(index).context("unsupported face"))
            .collect::<Result<Vec<_>>>()
    };
    let a = terms(a)?;
    let b = terms(b)?;
    let delta = |rank| {
        accurate_sum(
            a.iter()
                .map(|t| t.0[rank])
                .chain(b.iter().map(|t| -t.0[rank])),
        )
    };
    let volume = |terms: &[([f64; 2], [f64; 2])]| accurate_sum(terms.iter().map(|t| t.0[1]));
    let drift = |terms: &[([f64; 2], [f64; 2])]| {
        accurate_sum(terms.iter().flat_map(|t| [t.0[1], -t.0[0]])).map(f64::abs)
    };
    let error = |rank| accurate_sum(a.iter().chain(&b).map(|t| t.1[rank]));
    Ok(ComparisonSample {
        reference,
        volumes: [volume(&a)?, volume(&b)?],
        delta: delta(1)?,
        coarse_delta: delta(0)?,
        quadrature_drift: [drift(&a)?, drift(&b)?],
        estimated_quadrature_error: error(0)?.max(error(1)?),
    })
}

fn reference_spread(faces: &[Face]) -> Result<f64> {
    let mut max = 0.0_f64;
    for i in 0..4 {
        for j in 0..i {
            let terms = faces
                .iter()
                .map(|face| {
                    let a = face.flux.terms(i).context("unsupported face")?.0[1];
                    let b = face.flux.terms(j).context("unsupported face")?.0[1];
                    Ok([a, -b])
                })
                .collect::<Result<Vec<_>>>()?;
            max = max.max(accurate_sum(terms.into_iter().flatten())?.abs());
        }
    }
    Ok(max)
}

#[derive(Debug, Serialize)]
struct Measurement {
    samples: Vec<ComparisonSample>,
    reference_spread: [f64; 2],
    max_absolute_delta: f64,
    max_quadrature_drift: f64,
    max_estimated_quadrature_error: f64,
    observed_within_tolerance: bool,
}

#[derive(Debug, Serialize)]
pub(crate) struct SpatialVolumeReport {
    method: &'static str,
    absolute_tolerance_mm3: f64,
    boundary_tolerance_mm: f64,
    brep_valid: [bool; 2],
    unsupported_faces: [usize; 2],
    measurement: Option<Measurement>,
    a: Vec<Face>,
    b: Vec<Face>,
}

pub(crate) fn compare(
    kernel: &mut OcctKernel,
    a: ShapeHandle,
    b: ShapeHandle,
    absolute_tolerance: f64,
    boundary_tolerance: f64,
) -> Result<SpatialVolumeReport> {
    ensure!(
        [absolute_tolerance, boundary_tolerance]
            .iter()
            .all(|t| t.is_finite() && *t > 0.0),
        "invalid spatial tolerances"
    );
    let references = references(
        kernel.get_bounding_box_loose(a, true)?,
        kernel.get_bounding_box_loose(b, true)?,
    )?;
    let valid = [kernel.is_valid(a)?, kernel.is_valid(b)?];
    let a = faces(kernel, a, boundary_tolerance, references)?;
    let b = faces(kernel, b, boundary_tolerance, references)?;
    let unsupported = |faces: &[Face]| {
        faces
            .iter()
            .filter(|face| matches!(face.flux, Flux::Unsupported { .. }))
            .count()
    };
    let unsupported_faces = [unsupported(&a), unsupported(&b)];
    let measurement = if unsupported_faces == [0, 0] {
        let samples = references
            .into_iter()
            .enumerate()
            .map(|(i, r)| compare_sample(&a, &b, i, r))
            .collect::<Result<Vec<_>>>()?;
        let reference_spread = [reference_spread(&a)?, reference_spread(&b)?];
        let max_absolute_delta = samples
            .iter()
            .map(|s| s.delta.abs().max(s.coarse_delta.abs()))
            .fold(0.0, f64::max);
        let max_quadrature_drift = samples
            .iter()
            .map(|s| s.quadrature_drift[0].max(s.quadrature_drift[1]))
            .fold(0.0, f64::max);
        let max_estimated_quadrature_error = samples
            .iter()
            .map(|s| s.estimated_quadrature_error)
            .fold(0.0, f64::max);
        // Empirical integral agreement only; incidence/closure geometry
        // bounds are reported separately and are not interval volume bounds.
        let observed_within_tolerance = valid == [true, true]
            && samples.iter().all(|s| s.volumes.iter().all(|v| *v > 0.0))
            && max_absolute_delta
                .max(max_quadrature_drift)
                .max(reference_spread[0])
                .max(reference_spread[1])
                + max_estimated_quadrature_error
                <= absolute_tolerance;
        Some(Measurement {
            samples,
            reference_spread,
            max_absolute_delta,
            max_quadrature_drift,
            max_estimated_quadrature_error,
            observed_within_tolerance,
        })
    } else {
        None
    };
    Ok(SpatialVolumeReport {
        method: "independent elementary spatial boundary flux; diagnostic only",
        absolute_tolerance_mm3: absolute_tolerance,
        boundary_tolerance_mm: boundary_tolerance,
        brep_valid: valid,
        unsupported_faces,
        measurement,
        a,
        b,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_boxes_and_cylinders_give_correct_signed_spatial_volume() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let cube = kernel.make_box(2.0, 3.0, 4.0)?;
        let cylinder = kernel.make_cylinder(2.0, 3.0)?;
        for (shape, expected) in [(cube, 24.0), (cylinder, 12.0 * std::f64::consts::PI)] {
            let report = compare(&mut kernel, shape, shape, 1e-10, 2e-9)?;
            assert_eq!(report.unsupported_faces, [0, 0], "{report:?}");
            let measured = report.measurement.context("missing spatial volume")?;
            assert!(measured.observed_within_tolerance, "{measured:?}");
            for sample in measured.samples {
                assert!((sample.volumes[0] - expected).abs() < 1e-10, "{sample:?}");
                assert!(sample.delta.abs() < 1e-13);
            }
        }
        Ok(())
    }

    #[test]
    fn spatial_volume_handles_holes_partitions_and_small_real_changes() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let whole = kernel.make_box(2.0, 3.0, 4.0)?;
        let half = kernel.make_box(1.0, 3.0, 4.0)?;
        let right = kernel.translate(half, 1.0, 0.0, 0.0)?;
        let split = kernel.fuse(half, right)?;
        let report = compare(&mut kernel, whole, split, 1e-10, 2e-9)?;
        assert_eq!(report.unsupported_faces, [0, 0], "{report:?}");
        assert!(
            report
                .measurement
                .context("missing split volume")?
                .observed_within_tolerance
        );
        let altered = kernel.make_box(2.0, 3.0, 4.0 + 1e-8)?;
        let report = compare(&mut kernel, whole, altered, 1e-10, 2e-9)?;
        let measured = report.measurement.context("missing changed volume")?;
        assert!(!measured.observed_within_tolerance);
        for sample in measured.samples {
            assert!((sample.delta + 6e-8).abs() < 1e-10, "{sample:?}");
        }
        let outer = kernel.make_cylinder(2.0, 3.0)?;
        let inner = kernel.make_cylinder(1.0, 5.0)?;
        let inner = kernel.translate(inner, 0.0, 0.0, -1.0)?;
        let hollow = kernel.cut(outer, inner)?;
        let translated = kernel.translate(hollow, 67.75, 5.1, -9.0)?;
        let report = compare(&mut kernel, translated, translated, 1e-10, 2e-9)?;
        assert_eq!(report.unsupported_faces, [0, 0], "{report:?}");
        let measured = report.measurement.context("missing annular volume")?;
        assert!(measured.observed_within_tolerance, "{measured:?}");
        for sample in measured.samples {
            assert!(
                (sample.volumes[0] - 9.0 * std::f64::consts::PI).abs() < 1e-10,
                "{sample:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn closed_cylinder_with_reversed_occurrence_has_negative_flux() -> Result<()> {
        let carrier = Carrier::Cylinder {
            axis: [0., 0., 1.],
            origin: [0.; 3],
            radius: 2.,
        };
        let circle = |z| {
            spatial_boundary::circle_spans(
                [0., 0., z],
                1.,
                [0., 0., 1.],
                [1., 0., 0.],
                0.,
                std::f64::consts::TAU,
            )
        };
        let mut top = circle(1.5)?;
        let bottom = circle(0.)?;
        // Outward surface, theta-z chart: lower edge increases theta;
        // upper decreases theta. Both directed circles together close.
        top.reverse();
        for span in &mut top {
            span.reverse();
        }
        let mut spans = [top, bottom].concat();
        let positive = carrier.sample(&spans, 1., 1e-13, 0., [[0.; 3]; 4])?;
        assert!((positive.volume[0] - 8. * std::f64::consts::PI).abs() < 1e-10);
        for span in &mut spans {
            span.reverse();
        }
        let negative = carrier.sample(&spans, -1., 1e-13, 0., [[0.; 3]; 4])?;
        assert!((positive.volume[0] + negative.volume[0]).abs() < 1e-12);
        Ok(())
    }
}
