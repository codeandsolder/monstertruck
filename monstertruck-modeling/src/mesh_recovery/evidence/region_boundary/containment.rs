//! Analytic restrictions of a shared trim to independently fitted surface
//! carriers. A small vertex residual alone cannot establish containment.

use super::super::super::sample_fit::normal_basis;
use super::super::{CircularMeridianEvidence, MeshEvidence, PlanarSupport, SpatialTrackEvidence};
mod circular;
use super::{ContactCircleEvidence, ContactCurveEvidence, ContactLineEvidence};
use crate::{InnerSpace, Point3, Vector3};
use circular::CircularProfile;

/// Infinite surface carrier hypothesis, independently fitted to source faces.
///
/// Linear and circular meridians represent different rotational carriers.
#[derive(Clone, Copy, Debug)]
pub enum SurfaceCarrier {
    /// Oriented plane.
    Plane {
        /// Source-relative plane reference point.
        origin: Point3,
        /// Unit plane normal.
        normal: Vector3,
    },
    /// Surface of revolution with affine radius in axial coordinate.
    LinearRevolution {
        /// Reference point on the axis.
        origin: Point3,
        /// Unit axial direction.
        axis: Vector3,
        /// Profile radius at the reference point (may extrapolate below zero).
        radius_at_origin: f64,
        /// Change in radius per unit axial distance.
        slope: f64,
    },
    /// Surface of revolution whose meridian is circular (including spheres).
    CircularRevolution {
        /// Axis point at the meridian circle center's axial coordinate.
        axial_center: Point3,
        /// Unit axial direction.
        axis: Vector3,
        /// Signed radial coordinate of the meridian circle center.
        radial_center: f64,
        /// Radius of the meridian circle.
        profile_radius: f64,
    },
}

/// Numerical continuous-carrier residual, in mesh length units.
/// The bound is analytic (not a finite sample estimate), evaluated in f64;
/// it is not a directed-rounding interval arithmetic certificate.
#[derive(Clone, Copy, Debug)]
pub struct CarrierResidual {
    /// Upper bound on maximum absolute plane or meridian residual.
    pub upper_bound: f64,
    /// Maximum directly evaluated residual, a lower witness for the maximum.
    pub witnessed_max: f64,
    /// Minimum predicted radius anywhere along the tested carrier domain.
    /// Negative means the proposed trim extends beyond the selected profile
    /// sheet. Only defined for linear meridians.
    pub minimum_profile_radius: Option<f64>,
}

/// Independent carrier constraints for the competing shared trim curves.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceCurveAgreement {
    /// Residual of the finite endpoint line segment.
    pub line: Option<CarrierResidual>,
    /// Residual over the entire circular carrier, not just the observed arc.
    /// Exact containment of a nontrivial circular arc by an infinite plane
    /// or polynomial surface implies containment of the whole circle on its
    /// algebraic carrier; the chosen radial sheet is checked separately.
    pub circle: Option<CarrierResidual>,
}

impl SurfaceCarrier {
    // Signed implicit residual and gradient, in source coordinates.
    // Linear profile residual is normalized to unit gradient.
    pub(super) fn signed_constraint(self, point: Point3) -> Option<(Vector3, f64)> {
        match self {
            Self::Plane { origin, normal } => Some((normal, (point - origin).dot(normal))),
            Self::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } => {
                let delta = point - origin;
                let z = delta.dot(axis);
                let radial = delta - axis * z;
                let rho = radial.magnitude();
                if rho <= 0.0 {
                    return None;
                }
                let factor = 1.0_f64.hypot(slope);
                Some((
                    (radial / rho - axis * slope) / factor,
                    (rho - slope.mul_add(z, radius_at_origin)) / factor,
                ))
            }
            Self::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                profile_radius,
            } => {
                let delta = point - axial_center;
                if radial_center == 0.0 {
                    let rho = delta.magnitude();
                    return (rho > 0.0).then(|| (delta / rho, rho - profile_radius));
                }
                let z = delta.dot(axis);
                let radial = delta - axis * z;
                let rho = radial.magnitude();
                let distance = z.hypot(rho - radial_center);
                if rho <= 0.0 || distance <= 0.0 {
                    return None;
                }
                Some((
                    (axis * z + radial * ((rho - radial_center) / rho)) / distance,
                    distance - profile_radius,
                ))
            }
        }
    }

    /// Obtain the independent plane hypothesis of an original local support.
    #[must_use = "the plane is a candidate, not a recovered face"]
    pub const fn from_planar_support(support: &PlanarSupport) -> Self {
        Self::Plane {
            origin: support.origin,
            normal: support.normal,
        }
    }

    /// Use the independently fitted linear vertex-radius profile, rather than
    /// the apothem of the original polygonal support planes. An unresolved
    /// slope does not supply a linear-revolution hypothesis.
    #[must_use = "an unresolved radial profile has no carrier"]
    pub fn from_spatial_fit(fit: &SpatialTrackEvidence) -> Option<Self> {
        Some(Self::LinearRevolution {
            origin: fit.axis_origin,
            axis: fit.axis,
            radius_at_origin: fit.fitted_sample_radius_at_origin,
            slope: fit.fitted_sample_radius_slope?,
        })
    }

    /// Obtain the circular meridian hypothesis without discarding its carrier.
    #[must_use = "the meridian requires independent region validation"]
    pub const fn from_circular_meridian(fit: &CircularMeridianEvidence) -> Self {
        Self::CircularRevolution {
            axial_center: fit.axial_center,
            axis: fit.axis,
            radial_center: fit.radial_center,
            profile_radius: fit.profile_radius,
        }
    }

    /// Test the shared trim geometry against this independent surface.
    #[must_use = "inspect the residuals and radial-sheet evidence"]
    pub fn measure(self, curves: &ContactCurveEvidence) -> SurfaceCurveAgreement {
        SurfaceCurveAgreement {
            line: curves
                .line
                .as_ref()
                .and_then(|line| self.line_residual(line)),
            circle: curves
                .circle
                .as_ref()
                .and_then(|circle| self.circle_residual(circle)),
        }
    }

    fn line_residual(self, line: &ContactLineEvidence) -> Option<CarrierResidual> {
        let [start, end] = line.endpoints;
        self.measure_segment(start, end)
    }

    /// Bound the continuous residual on a finite source chord.
    /// This does not require the chord to be an intended analytic trim.
    #[must_use = "the witnessed and upper bounds have different meanings"]
    pub fn measure_segment(self, start: Point3, end: Point3) -> Option<CarrierResidual> {
        match self {
            Self::Plane { origin, normal } if unit(normal) => {
                let error = (start - origin)
                    .dot(normal)
                    .abs()
                    .max((end - origin).dot(normal).abs());
                residual(error, error, None)
            }
            Self::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } if unit(axis) && radius_at_origin.is_finite() && slope.is_finite() => {
                let a = RadialTerm::new(start - origin, axis, radius_at_origin, slope);
                let b = RadialTerm::new(end - start, axis, 0.0, slope);
                let c0 = a.form(a);
                let c1 = 2.0 * a.form(b);
                let c2 = b.form(b);
                let mut maximum = c0.abs().max((c0 + c1 + c2).abs());
                let minimum_radius = a.expected.min(a.expected + b.expected);
                let mut witness = a.error().max(a.add_scaled(b, 1.0).error());
                if c2 != 0.0 {
                    let extremum = -c1 / (2.0 * c2);
                    if (0.0..=1.0).contains(&extremum) {
                        maximum = maximum.max(c2.mul_add(extremum, c1).mul_add(extremum, c0).abs());
                        witness = witness.max(a.add_scaled(b, extremum).error());
                    }
                }
                // The minimum radial distance may occur inside a chord.
                let denominator = b.radial.magnitude2();
                if denominator > 0.0 {
                    let t = (-a.radial.dot(b.radial) / denominator).clamp(0.0, 1.0);
                    witness = witness.max(a.add_scaled(b, t).error());
                }
                radial_residual(maximum, witness, minimum_radius)
            }
            Self::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                profile_radius,
            } => CircularProfile::new(axial_center, axis, radial_center, profile_radius)?
                .line(start, end),
            _ => None,
        }
    }

    /// Measure all finite polygon-edge chords in the selected original
    /// supports, including internal tessellation edges. The maximum witness
    /// is a lower bound, and the maximum upper bound covers every chord.
    /// This is not a bound over polygon interiors or the reverse direction.
    /// When a continuous bound straddles the specified chord budget,
    /// subdivide the chord rather than confuse a loose bound with an error.
    #[must_use = "chord approximation is distinct from vertex interpolation"]
    pub fn measure_source_chords(
        self,
        evidence: &MeshEvidence,
        supports: &[usize],
        chord_budget: f64,
    ) -> Option<CarrierResidual> {
        if chord_budget.is_nan() || chord_budget < 0.0 {
            return None;
        }
        let mut result = None::<CarrierResidual>;
        for &id in supports {
            for &polygon in &evidence.planar_supports.get(id)?.source_face_indices {
                let vertices = evidence.source_polygon_vertices(polygon)?;
                for (&a, &b) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
                    let measured = self.measure_segment_to_budget(
                        *evidence.positions.get(a)?,
                        *evidence.positions.get(b)?,
                        chord_budget,
                        8,
                    )?;
                    if let Some(previous) = &mut result {
                        previous.upper_bound = previous.upper_bound.max(measured.upper_bound);
                        previous.witnessed_max = previous.witnessed_max.max(measured.witnessed_max);
                        if let Some(radius) = measured.minimum_profile_radius {
                            previous.minimum_profile_radius = Some(
                                previous
                                    .minimum_profile_radius
                                    .map_or(radius, |old| old.min(radius)),
                            );
                        }
                    } else {
                        result = Some(measured);
                    }
                }
            }
        }
        result
    }

    /// Measure the filled original source facets of a circular revolution.
    /// A perimeter-only chord certificate does not bound polygon interiors.
    /// All original triangles (including the mesh's polygon fan
    /// triangulation) are checked against the same absolute budget.
    #[must_use = "the full-facet bound is distinct from source chord evidence"]
    pub fn measure_circular_source_facets(
        self,
        evidence: &MeshEvidence,
        supports: &[usize],
        budget: f64,
    ) -> Option<CarrierResidual> {
        let Self::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } = self
        else {
            return None;
        };
        let profile = CircularProfile::new(axial_center, axis, radial_center, profile_radius)?;
        let mut result = None::<CarrierResidual>;
        for &id in supports {
            for &polygon in &evidence.planar_supports.get(id)?.source_face_indices {
                let vertices = evidence.source_polygon_vertices(polygon)?;
                let points = vertices
                    .iter()
                    .map(|&id| evidence.positions.get(id).copied())
                    .collect::<Option<Vec<_>>>()?;
                let measured = profile.facet(&points, budget)?;
                result = Some(result.map_or(measured, |previous| CarrierResidual {
                    upper_bound: previous.upper_bound.max(measured.upper_bound),
                    witnessed_max: previous.witnessed_max.max(measured.witnessed_max),
                    minimum_profile_radius: None,
                }));
            }
        }
        result
    }

    fn measure_segment_to_budget(
        self,
        start: Point3,
        end: Point3,
        budget: f64,
        depth: usize,
    ) -> Option<CarrierResidual> {
        let measured = self.measure_segment(start, end)?;
        if depth == 0 || measured.upper_bound <= budget || measured.witnessed_max > budget {
            return Some(measured);
        }
        let middle = start + (end - start) * 0.5;
        let left = self.measure_segment_to_budget(start, middle, budget, depth - 1)?;
        let right = self.measure_segment_to_budget(middle, end, budget, depth - 1)?;
        Some(CarrierResidual {
            upper_bound: left.upper_bound.max(right.upper_bound),
            witnessed_max: left.witnessed_max.max(right.witnessed_max),
            minimum_profile_radius: left
                .minimum_profile_radius
                .zip(right.minimum_profile_radius)
                .map(|(a, b)| a.min(b)),
        })
    }

    /// Finite observed circle arc on the selected radial sheet.
    /// Keep whole-circle carrier comparisons separate from this trim proof.
    pub(super) fn measure_circle_arc(
        self,
        circle: &ContactCircleEvidence,
        first: Point3,
    ) -> Option<CarrierResidual> {
        let Self::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } = self
        else {
            return self.circle_residual(circle);
        };
        if !unit(circle.normal) || !circle.radius.is_finite() || circle.radius <= 0.0 {
            return None;
        }
        let relative = first - circle.center;
        let radial = relative - circle.normal * relative.dot(circle.normal);
        let length = radial.magnitude();
        if length <= 0.0 || !length.is_finite() {
            return None;
        }
        let u = radial * (circle.radius / length);
        let v = circle.normal.cross(u);
        CircularProfile::new(axial_center, axis, radial_center, profile_radius)?.circle_over(
            circle.center,
            u,
            v,
            circle.radius,
            circle.signed_sweep_radians,
        )
    }

    fn circle_residual(self, circle: &ContactCircleEvidence) -> Option<CarrierResidual> {
        if !unit(circle.normal) || !circle.radius.is_finite() || circle.radius <= 0.0 {
            return None;
        }
        let (basis_u, basis_v) = normal_basis(circle.normal);
        let (basis_u, basis_v) = (basis_u * circle.radius, basis_v * circle.radius);
        match self {
            Self::Plane { origin, normal } if unit(normal) => {
                let error = (circle.center - origin).dot(normal).abs()
                    + basis_u.dot(normal).hypot(basis_v.dot(normal));
                residual(error, error, None)
            }
            Self::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } if unit(axis) && radius_at_origin.is_finite() && slope.is_finite() => {
                let a = RadialTerm::new(circle.center - origin, axis, radius_at_origin, slope);
                let b = RadialTerm::new(basis_u, axis, 0.0, slope);
                let c = RadialTerm::new(basis_v, axis, 0.0, slope);
                // Restrict Q = radial² - expected_radius² to
                // a + b cos(t) + c sin(t). Its Fourier series has degree two.
                let bb = b.form(b);
                let cc = c.form(c);
                let constant = a.form(a) + bb.midpoint(cc);
                let first = 2.0 * a.form(b).hypot(a.form(c));
                let second = ((bb - cc) * 0.5).hypot(b.form(c));
                let maximum = constant.abs() + first + second;
                let minimum_radius = a.expected - b.expected.hypot(c.expected);
                let witness = (0..8)
                    .map(|i| {
                        let (sin, cos) = (f64::from(i) * std::f64::consts::FRAC_PI_4).sin_cos();
                        a.add_scaled(b, cos).add_scaled(c, sin).error()
                    })
                    .fold(0.0_f64, f64::max);
                radial_residual(maximum, witness, minimum_radius)
            }
            Self::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                profile_radius,
            } => CircularProfile::new(axial_center, axis, radial_center, profile_radius)?.circle(
                circle.center,
                basis_u,
                basis_v,
                circle.radius,
            ),
            _ => None,
        }
    }
}

pub(super) fn unit(axis: Vector3) -> bool {
    (0.0..=64.0 * f64::EPSILON).contains(&(axis.magnitude2() - 1.0).abs())
}

fn residual(
    bound: f64,
    witness: f64,
    minimum_profile_radius: Option<f64>,
) -> Option<CarrierResidual> {
    if !bound.is_finite()
        || !witness.is_finite()
        || minimum_profile_radius.is_some_and(|r| !r.is_finite())
    {
        return None;
    }
    Some(CarrierResidual {
        upper_bound: bound.max(witness),
        witnessed_max: witness,
        minimum_profile_radius,
    })
}

/// Difference-of-squares bound: |r-e| = |r²-e²|/(r+e) when e>0.
/// Near an apex use |r-|e|| ≤ sqrt(|r²-e²|); for a negative expected
/// radius also account for 2|e| and report the sheet crossing separately.
fn radial_residual(
    squared_bound: f64,
    witness: f64,
    minimum_radius: f64,
) -> Option<CarrierResidual> {
    let root = squared_bound.sqrt();
    let bound = if minimum_radius > 0.0 {
        root.min(squared_bound / minimum_radius)
    } else {
        2.0f64.mul_add(-minimum_radius, root)
    };
    residual(bound, witness, Some(minimum_radius))
}

#[derive(Clone, Copy)]
struct RadialTerm {
    radial: Vector3,
    expected: f64,
}

impl RadialTerm {
    fn new(delta: Vector3, axis: Vector3, radius: f64, slope: f64) -> Self {
        let axial = delta.dot(axis);
        Self {
            radial: delta - axis * axial,
            expected: slope.mul_add(axial, radius),
        }
    }

    fn form(self, other: Self) -> f64 {
        self.expected
            .mul_add(-other.expected, self.radial.dot(other.radial))
    }

    fn add_scaled(self, other: Self, scale: f64) -> Self {
        Self {
            radial: self.radial + other.radial * scale,
            expected: other.expected.mul_add(scale, self.expected),
        }
    }

    fn error(self) -> f64 { (self.radial.magnitude() - self.expected).abs() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    fn line(start: Point3, end: Point3) -> ContactLineEvidence {
        ContactLineEvidence {
            endpoints: [start, end],
            direction: (end - start).normalize(),
            chord_length: (end - start).magnitude(),
            max_vertex_offset: 0.0,
            max_reverse_step: 0.0,
            exactly_collinear: true,
        }
    }

    fn circle(center: Point3, normal: Vector3, radius: f64) -> ContactCircleEvidence {
        ContactCircleEvidence {
            center,
            normal,
            radius,
            plane_span_strength: 1.0,
            max_plane_residual: 0.0,
            max_vertex_radial_residual: 0.0,
            max_chord_radial_residual: 0.0,
            signed_sweep_radians: TAU,
            unsigned_sweep_radians: TAU,
            max_step_radians: 0.1,
            zero_angular_steps: 0,
        }
    }

    #[test]
    fn independent_plane_and_cone_contain_same_circle() -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(1.0e6, -2.0e6, 3.0);
        let axis = Vector3::unit_z();
        let center = origin + axis * 4.0;
        let profile = SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: 2.0,
            slope: 0.5,
        };
        let plane = SurfaceCarrier::Plane {
            origin: center,
            normal: axis,
        };
        let curves = ContactCurveEvidence {
            source_vertex_count: 16,
            closed: true,
            source_length: TAU * 4.0,
            line: None,
            circle: Some(circle(center, axis, 4.0)),
        };
        for carrier in [plane, profile] {
            let measured = carrier
                .measure(&curves)
                .circle
                .ok_or("missing circle constraint")?;
            assert!(measured.upper_bound < 1.0e-12);
        }
        let shifted = circle(center + Vector3::unit_x() * 0.125, axis, 4.0);
        let residual = profile
            .circle_residual(&shifted)
            .ok_or("missing circle residual")?;
        assert!(residual.witnessed_max >= 0.1249);
        assert!(residual.upper_bound >= residual.witnessed_max);
        Ok(())
    }

    #[test]
    fn endpoint_containment_does_not_prove_a_cylindrical_chord()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let carrier = SurfaceCarrier::LinearRevolution {
            origin,
            axis: Vector3::unit_z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let chord = line(Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let residual = carrier
            .line_residual(&chord)
            .ok_or("missing line residual")?;
        assert!((residual.witnessed_max - 1.0).abs() < 1.0e-12);
        assert!((residual.upper_bound - 1.0).abs() < 1.0e-12);
        let generator = line(Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 5.0));
        assert!(
            carrier
                .line_residual(&generator)
                .ok_or("missing generator")?
                .upper_bound
                < 1.0e-12
        );
        Ok(())
    }

    #[test]
    fn source_chord_bounds_distinguish_planar_facets_from_exact_vertex_interpolants()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{EvidenceOptions, analyze_mesh_evidence};
        use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
        let positions = (0..3)
            .map(|i| {
                let (sin, cos) = (TAU * f64::from(i) / 3.0).sin_cos();
                Point3::new(21.0 * cos, 21.0 * sin, 0.0)
            })
            .collect::<Vec<_>>();
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: 20.0,
            profile_radius: 1.0,
        };
        for &point in &evidence.positions {
            assert!(
                carrier
                    .measure_segment(point, point)
                    .ok_or("vertex residual")?
                    .witnessed_max
                    < 1.0e-12
            );
        }
        let residual = carrier
            .measure_source_chords(&evidence, &[0], f64::INFINITY)
            .ok_or("chord residual")?;
        assert!(residual.witnessed_max > 8.0);
        assert!(residual.upper_bound >= residual.witnessed_max);
        let plane = SurfaceCarrier::from_planar_support(&evidence.planar_supports[0]);
        assert!(
            plane
                .measure_source_chords(&evidence, &[0], f64::INFINITY)
                .ok_or("plane residual")?
                .upper_bound
                < 1.0e-12
        );
        Ok(())
    }

    #[test]
    fn squared_cone_equation_alone_cannot_select_the_radial_sheet()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let cone = SurfaceCarrier::LinearRevolution {
            origin,
            axis: Vector3::unit_z(),
            radius_at_origin: 0.0,
            slope: 1.0,
        };
        let wrong_sheet = circle(Point3::new(0.0, 0.0, -2.0), Vector3::unit_z(), 2.0);
        let residual = cone
            .circle_residual(&wrong_sheet)
            .ok_or("missing sheet evidence")?;
        assert!(residual.minimum_profile_radius.is_some_and(|r| r < 0.0));
        assert!((residual.witnessed_max - 4.0).abs() < 1.0e-12);
        assert!((residual.upper_bound - 4.0).abs() < 1.0e-12);
        let crossing = line(Point3::new(-2.0, 0.0, -2.0), Point3::new(2.0, 0.0, 2.0));
        let residual = cone
            .line_residual(&crossing)
            .ok_or("missing sheet evidence")?;
        assert!(residual.minimum_profile_radius.is_some_and(|r| r < 0.0));
        assert!(residual.upper_bound >= 4.0);
        Ok(())
    }

    #[test]
    fn full_circle_harmonic_bound_covers_off_axis_and_tilted_cases()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(2.0, -3.0, 1.0);
        let axis = Vector3::new(0.1, 0.2, 1.0).normalize();
        let normal = Vector3::new(-0.2, 0.5, 0.8).normalize();
        let candidate = circle(Point3::new(1.5, -2.8, 4.0), normal, 1.3);
        let carrier = SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: 1.0,
            slope: 0.25,
        };
        let measured = carrier
            .circle_residual(&candidate)
            .ok_or("missing circle residual")?;
        let (u, v) = normal_basis(normal);
        for i in 0..4096 {
            let (sin, cos) = (TAU * f64::from(i) / 4096.0).sin_cos();
            let delta = candidate.center - origin + (u * cos + v * sin) * candidate.radius;
            let sample = RadialTerm::new(delta, axis, 1.0, 0.25);
            assert!(sample.error() <= measured.upper_bound + 1.0e-12);
            assert!(
                measured
                    .minimum_profile_radius
                    .is_some_and(|r| r <= sample.expected + 1.0e-12)
            );
        }
        Ok(())
    }
}
