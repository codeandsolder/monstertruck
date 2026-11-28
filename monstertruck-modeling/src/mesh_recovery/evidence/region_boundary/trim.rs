//! Shared analytic trim hypotheses over contiguous original source paths.
//! Segmentation never invents a second copy of the same edge on either side.

use super::curve_evidence::{measure_path_curves, measure_placed_circle};
use super::{
    ContactCircleEvidence, ContactConicEvidence, ContactLineEvidence, ContactLinearSectionEvidence,
    ContactToroidalSectionEvidence, RegionContactPath, SurfaceCarrier,
};
use crate::Point3;
use std::f64::consts::TAU;

/// Local fit budget and bounded search configuration.
#[derive(Clone, Copy, Debug)]
pub struct SharedTrimOptions {
    /// Absolute source sample/continuous carrier error budget.
    pub absolute_tolerance: f64,
    /// Budget relative to observed source path length and carrier radius.
    pub relative_tolerance: f64,
    /// Longest interval examined during mixed-path segmentation. Whole
    /// original paths are always tested first, regardless of this limit.
    pub max_interval_edges: usize,
}

/// Analytic curve hypotheses certified against the supplied surface carriers.
#[derive(Clone, Copy, Debug)]
pub enum SharedTrimGeometry {
    /// Finite endpoint line with independent source residual measurements.
    Line(ContactLineEvidence),
    /// Monotone circular arc or full ring with independent sample evidence.
    Circle(ContactCircleEvidence),
    /// Finite rational conic branch with independently certified carrier containment.
    Conic(ContactConicEvidence),
    /// Regular plane/ring-torus section with continuously certified native approximation.
    ToroidalSection(ContactToroidalSectionEvidence),
    /// Regular shared nonconic section of two affine-radius revolutions.
    LinearSection(ContactLinearSectionEvidence),
    /// No accepted analytic hypothesis for these original edges.
    Unresolved,
}

#[derive(Clone, Copy)]
pub(super) enum SectionBranch {
    Toroidal(ContactToroidalSectionEvidence),
    Linear(ContactLinearSectionEvidence),
}

impl SectionBranch {
    pub(super) fn projected(self, point: Point3) -> Option<Point3> {
        match self {
            Self::Toroidal(section) => section.projected(point),
            Self::Linear(section) => section.projected(point),
        }
    }
}

impl SharedTrimGeometry {
    pub(super) const fn section(self) -> Option<SectionBranch> {
        match self {
            Self::ToroidalSection(section) => Some(SectionBranch::Toroidal(section)),
            Self::LinearSection(section) => Some(SectionBranch::Linear(section)),
            _ => None,
        }
    }
}

/// One contiguous piece of a globally shared directed source path.
#[derive(Clone, Copy, Debug)]
pub struct SharedTrimPiece {
    /// First original edge index in the input path (cyclic for closed paths).
    pub first_edge: usize,
    /// Number of successive original edges, preserving original winding.
    pub edge_count: usize,
    /// One common geometry for both surface sides.
    pub geometry: SharedTrimGeometry,
}

/// Shared path decomposition. Unresolved edges remain represented.
#[derive(Debug)]
pub struct SharedTrimSegmentation {
    /// Ordered disjoint contiguous original edge ranges.
    pub pieces: Vec<SharedTrimPiece>,
    /// Whether the original input path closes on its first source vertex.
    pub closed: bool,
    /// True when mixed-path search was restricted by the interval budget.
    pub limited: bool,
}

impl RegionContactPath {
    /// Fit and segment the shared source path against every adjacent carrier.
    /// A sparse circle derived from the surfaces or a three-point fit requires distinct incident carrier
    /// constraints. Otherwise at least four samples and a 60-degree
    /// sampling envelope are required. All circular candidates are tested
    /// over their whole carrier, not just their sampled vertices. Missing adjacent carriers are not
    /// silently ignored. Mixed paths are optimized for original edge
    /// coverage first, then for the fewest pieces.
    #[must_use = "unresolved source edges must remain represented"]
    pub fn segment_analytic_trims(
        &self,
        points: &[Point3],
        carriers: &[SurfaceCarrier],
        options: SharedTrimOptions,
    ) -> Option<SharedTrimSegmentation> {
        if carriers.len() != 1 + usize::from(self.regions.1.is_some())
            || !options.absolute_tolerance.is_finite()
            || options.absolute_tolerance <= 0.0
            || !options.relative_tolerance.is_finite()
            || options.relative_tolerance <= 0.0
            || options.max_interval_edges == 0
        {
            return None;
        }
        let ids = self.ordered_vertices()?;
        if ids.iter().any(|&id| id >= points.len()) {
            return None;
        }
        let count = self.edges.len();
        let closed = ids.first() == ids.last();
        if let Some(geometry) = accepted(points, &ids, carriers, options) {
            return Some(SharedTrimSegmentation {
                pieces: vec![SharedTrimPiece {
                    first_edge: 0,
                    edge_count: count,
                    geometry,
                }],
                closed,
                limited: false,
            });
        }
        let mut pieces = segment_ranges(points, &ids, carriers, options);
        // A canonical cyclic path may start in the middle of one analytic
        // arc or line. Test the wraparound union of both end pieces.
        if self.unanchored_cycle && closed && pieces.len() > 1 {
            merge_cyclic_ends(points, &ids, carriers, options, &mut pieces);
        }
        Some(SharedTrimSegmentation {
            pieces,
            closed,
            limited: count > options.max_interval_edges,
        })
    }
}

fn accepted(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    options: SharedTrimOptions,
) -> Option<SharedTrimGeometry> {
    let fit = measure_path_curves(points, ids)?;
    let tolerance = options
        .absolute_tolerance
        .min(options.relative_tolerance * fit.source_length);
    if let Some(line) = fit.line
        && line.max_vertex_offset <= tolerance
        && line.max_reverse_step <= tolerance
        && carriers.iter().all(|carrier| {
            carrier
                .measure_segment(line.endpoints[0], line.endpoints[1])
                .is_some_and(|r| r.upper_bound <= tolerance)
        })
    {
        return Some(SharedTrimGeometry::Line(line));
    }
    let first = points[ids[0]];
    let span = ids
        .iter()
        .map(|&id| (points[id] - first).norm())
        .fold(0.0_f64, f64::max);
    // A free planar curve has no independent intersection fixing its
    // carrier. Its endpoints must therefore be constraints of the fit,
    // not independently adjustable predictions at a shared junction.
    let circle = if carriers
        .iter()
        .all(|c| matches!(c, SurfaceCarrier::Plane { .. }))
    {
        fit.circle
            .and_then(|observed| endpoint_circle(points, ids, carriers[0], observed))
    } else {
        fit.circle
    };
    if let Some(circle) = circle
        && accepts_circle(
            circle,
            fit,
            carriers,
            tolerance.min(options.relative_tolerance * circle.radius),
            span,
            first,
            false,
        )
    {
        return Some(SharedTrimGeometry::Circle(circle));
    }

    // A two-point path cannot fit a circle, but the independently recovered
    // surface intersection may supply it. Reject competing distinct circles
    // rather than arbitrarily choose one of the surviving analytic families.
    let mut selected = None;
    for circle in super::trim_intersections::supplied_circles(points, ids, carriers, tolerance) {
        if !accepts_circle(
            circle,
            fit,
            carriers,
            tolerance.min(options.relative_tolerance * circle.radius),
            span,
            first,
            true,
        ) {
            continue;
        }
        if let Some(previous) = selected {
            if !super::trim_intersections::same_circle(&previous, &circle, tolerance) {
                return None;
            }
        } else {
            selected = Some(circle);
        }
    }
    selected
        .map(SharedTrimGeometry::Circle)
        .or_else(|| {
            super::conic::supplied_conic(points, ids, carriers, tolerance)
                .map(SharedTrimGeometry::Conic)
        })
        .or_else(|| {
            super::toroidal_section::supplied_section(points, ids, carriers, tolerance)
                .map(SharedTrimGeometry::ToroidalSection)
        })
        .or_else(|| {
            super::linear_section::supplied_section(points, ids, carriers, tolerance)
                .map(SharedTrimGeometry::LinearSection)
        })
}

// Algebraic least squares in the one-dimensional family of circles
// through both source endpoints. Unlike an unconstrained three-parameter
// circle fit this cannot create a geometric gap between adjacent pieces.
fn endpoint_circle(
    points: &[Point3],
    ids: &[usize],
    carrier: SurfaceCarrier,
    observed: ContactCircleEvidence,
) -> Option<ContactCircleEvidence> {
    let SurfaceCarrier::Plane { origin, normal } = carrier else {
        return None;
    };
    let project = |p: Point3| p - normal * (p - origin).dot(&normal);
    let first = project(points[*ids.first()?]);
    let last = project(points[*ids.last()?]);
    let (center, radius) = if ids.first() == ids.last() {
        let center = project(observed.center);
        (center, (first - center).norm())
    } else {
        let delta = last - first;
        let length = delta.norm();
        if length <= 0.0 {
            return None;
        }
        let mid = first + delta * 0.5;
        let direction = normal.cross(&(delta / length));
        let half_squared = 0.25 * length * length;
        let (numerator, denominator) = ids.iter().fold((0.0_f64, 0.0_f64), |(num, den), &id| {
            let local = project(points[id]) - mid;
            let side = 2.0 * local.dot(&direction);
            (
                side.mul_add(local.norm_squared() - half_squared, num),
                side.mul_add(side, den),
            )
        });
        if denominator <= f64::EPSILON * length * length {
            return None;
        }
        let offset = numerator / denominator;
        (mid + direction * offset, (length * 0.5).hypot(offset))
    };
    measure_placed_circle(
        points,
        ids,
        center,
        normal,
        radius,
        observed.plane_span_strength,
    )
}

fn accepts_circle(
    circle: ContactCircleEvidence,
    fit: super::ContactCurveEvidence,
    carriers: &[SurfaceCarrier],
    tolerance: f64,
    span: f64,
    first: Point3,
    supplied: bool,
) -> bool {
    let angular = tolerance / circle.radius;
    let independently_constrained = independently_constrained_circle(&circle, carriers, tolerance);
    // An open path still has one independent surface. Three noncollinear
    // source points can test a unique circle supplied by that surface,
    // but two points alone cannot choose a trim on an open face.
    let source_corroborates = supplied
        && carriers.len() == 1
        && fit.source_vertex_count >= 3
        && fit.circle.is_some_and(|observed| {
            observed.plane_span_strength * span > tolerance
                && super::trim_intersections::same_circle(&observed, &circle, tolerance)
        });
    let corroborated = independently_constrained || source_corroborates;
    let max_step = if corroborated {
        std::f64::consts::FRAC_PI_2
    } else {
        std::f64::consts::FRAC_PI_3
    };
    if supplied && !corroborated {
        return false;
    }
    if (!corroborated && fit.source_vertex_count < 4)
        || (!supplied && circle.plane_span_strength * span <= tolerance)
        || circle.max_step_radians > max_step + angular
        || circle.max_plane_residual > tolerance
        || circle.max_vertex_radial_residual > tolerance
        || circle.unsigned_sweep_radians - circle.signed_sweep_radians.abs() > angular
        || circle.unsigned_sweep_radians > TAU + angular
        || circle.zero_angular_steps > 0
    {
        return false;
    }
    carriers.iter().all(|carrier| {
        carrier
            .measure_circle_arc(&circle, first)
            .is_some_and(|r| r.upper_bound <= tolerance)
    })
}

// A distinct surface constraint must be witnessed away from the proposed
// common circle. Merely having two region IDs does not supply extra geometry.
// Compare absolute residuals, independent of normal direction and axis sign.
fn independently_constrained_circle(
    circle: &ContactCircleEvidence,
    carriers: &[SurfaceCarrier],
    tolerance: f64,
) -> bool {
    let [left, right] = carriers else {
        return false;
    };
    let (u, v) = super::super::super::sample_fit::normal_basis(circle.normal);
    let radius = circle.radius;
    let directions = [
        crate::Vector3::new(0.0, 0.0, 0.0),
        circle.normal,
        -circle.normal,
        u * 0.5,
        -u * 0.5,
        v * 0.5,
        -v * 0.5,
    ];
    directions.into_iter().any(|direction| {
        let point = circle.center + direction * radius;
        left.measure_segment(point, point)
            .zip(right.measure_segment(point, point))
            .is_some_and(|(a, b)| (a.witnessed_max - b.witnessed_max).abs() > tolerance * 2.0)
    })
}

#[derive(Clone, Copy)]
struct Step {
    unresolved: usize,
    pieces: usize,
    end: usize,
    geometry: SharedTrimGeometry,
}

fn segment_ranges(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    options: SharedTrimOptions,
) -> Vec<SharedTrimPiece> {
    let count = ids.len() - 1;
    let mut best = vec![
        Step {
            unresolved: 0,
            pieces: 0,
            end: count,
            geometry: SharedTrimGeometry::Unresolved
        };
        count + 1
    ];
    for start in (0..count).rev() {
        let tail = best[start + 1];
        let mut selected = Step {
            unresolved: tail.unresolved + 1,
            pieces: tail.pieces + 1,
            end: start + 1,
            geometry: SharedTrimGeometry::Unresolved,
        };
        for end in start + 1..=count.min(start.saturating_add(options.max_interval_edges)) {
            let tail = best[end];
            if (tail.unresolved, tail.pieces + 1) > (selected.unresolved, selected.pieces) {
                continue;
            }
            let Some(geometry) = accepted(points, &ids[start..=end], carriers, options) else {
                continue;
            };
            if (tail.unresolved, tail.pieces + 1) <= (selected.unresolved, selected.pieces) {
                selected = Step {
                    unresolved: tail.unresolved,
                    pieces: tail.pieces + 1,
                    end,
                    geometry,
                };
            }
        }
        best[start] = selected;
    }
    let mut pieces: Vec<SharedTrimPiece> = Vec::new();
    let mut start = 0;
    while start < count {
        let step = best[start];
        if let Some(previous) = pieces.last_mut()
            && matches!(previous.geometry, SharedTrimGeometry::Unresolved)
            && matches!(step.geometry, SharedTrimGeometry::Unresolved)
        {
            previous.edge_count += step.end - start;
        } else {
            pieces.push(SharedTrimPiece {
                first_edge: start,
                edge_count: step.end - start,
                geometry: step.geometry,
            });
        }
        start = step.end;
    }
    pieces
}

fn merge_cyclic_ends(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    options: SharedTrimOptions,
    pieces: &mut Vec<SharedTrimPiece>,
) {
    let (Some(first), Some(last)) = (pieces.first(), pieces.last()) else {
        return;
    };
    let count = ids.len() - 1;
    let length = first.edge_count + last.edge_count;
    let first_edge = last.first_edge;
    let wrapped = (0..=length)
        .map(|i| ids[(first_edge + i) % count])
        .collect::<Vec<_>>();
    let geometry = if matches!(first.geometry, SharedTrimGeometry::Unresolved)
        && matches!(last.geometry, SharedTrimGeometry::Unresolved)
    {
        SharedTrimGeometry::Unresolved
    } else if let Some(geometry) = accepted(points, &wrapped, carriers, options) {
        geometry
    } else {
        return;
    };
    pieces.remove(0);
    pieces.pop();
    pieces.push(SharedTrimPiece {
        first_edge,
        edge_count: length,
        geometry,
    });
}

#[cfg(test)]
mod tests {
    use super::super::RegionBoundaryEdge;
    use super::*;
    use crate::Vector3;

    fn path(count: usize, closed: bool) -> RegionContactPath {
        RegionContactPath {
            regions: (0, Some(1)),
            edges: (0..count)
                .map(|id| RegionBoundaryEdge {
                    vertices: [id, if closed { (id + 1) % count } else { id + 1 }],
                    interior_polygon: id,
                    exterior_polygon: Some(id + count),
                })
                .collect(),
            unanchored_cycle: closed,
        }
    }

    fn options() -> SharedTrimOptions {
        SharedTrimOptions {
            absolute_tolerance: 1.0e-8,
            relative_tolerance: 1.0e-4,
            max_interval_edges: 128,
        }
    }

    #[test]
    fn approximate_planar_arcs_share_exact_source_junctions()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = (0..9)
            .map(|i| {
                let (sin, cos) = (f64::from(i) * 0.1).sin_cos();
                let radius = if i == 0 || i == 4 || i == 8 {
                    1.0
                } else {
                    1.0 + 1.0e-6 * f64::from(i % 3)
                };
                Point3::new(radius * cos, radius * sin, 0.0)
            })
            .collect::<Vec<_>>();
        let plane = SurfaceCarrier::Plane {
            origin: points[0],
            normal: crate::Vector3::z(),
        };
        let options = SharedTrimOptions {
            absolute_tolerance: 1.0e-5,
            relative_tolerance: 1.0e-4,
            max_interval_edges: 128,
        };
        let mut circles = Vec::new();
        for ids in [vec![0, 1, 2, 3, 4], vec![4, 5, 6, 7, 8]] {
            let Some(SharedTrimGeometry::Circle(circle)) =
                accepted(&points, &ids, &[plane], options)
            else {
                return Err("circle".into());
            };
            for id in [ids[0], *ids.last().ok_or("end")?] {
                assert!(((points[id] - circle.center).norm() - circle.radius).abs() < 1.0e-12);
            }
            circles.push(circle);
        }
        assert!((circles[0].center - circles[1].center).norm() > 1.0e-8);
        Ok(())
    }

    #[test]
    fn shared_contact_segments_mixed_curves_and_checks_both_surfaces()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut points = vec![Point3::new(2.0, 0.0, 0.0)];
        points.extend((0..=8).map(|i| {
            let (sin, cos) = (std::f64::consts::FRAC_PI_2 * f64::from(i) / 8.0).sin_cos();
            Point3::new(cos, sin, 0.0)
        }));
        points.push(Point3::new(0.0, 2.0, 0.0));
        let path = path(points.len() - 1, false);
        let plane = SurfaceCarrier::Plane {
            origin: points[0],
            normal: Vector3::z(),
        };
        let result = path
            .segment_analytic_trims(&points, &[plane, plane], options())
            .ok_or("segmentation")?;
        assert_eq!(result.pieces.len(), 3);
        assert!(matches!(
            result.pieces[0].geometry,
            SharedTrimGeometry::Line(_)
        ));
        assert!(matches!(
            result.pieces[1].geometry,
            SharedTrimGeometry::Circle(_)
        ));
        assert_eq!(result.pieces[1].edge_count, 8);
        assert!(matches!(
            result.pieces[2].geometry,
            SharedTrimGeometry::Line(_)
        ));
        let cylinder = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let result = path
            .segment_analytic_trims(&points, &[plane, cylinder], options())
            .ok_or("segmentation")?;
        assert_eq!(result.pieces.len(), 3);
        assert!(matches!(
            result.pieces[0].geometry,
            SharedTrimGeometry::Unresolved
        ));
        assert!(matches!(
            result.pieces[1].geometry,
            SharedTrimGeometry::Circle(_)
        ));
        assert!(matches!(
            result.pieces[2].geometry,
            SharedTrimGeometry::Unresolved
        ));
        assert!(
            path.segment_analytic_trims(&points, &[plane], options())
                .is_none()
        );
        Ok(())
    }

    #[test]
    fn two_samples_use_independent_axial_circle() -> Result<(), Box<dyn std::error::Error>> {
        let points = [Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)];
        let plane = SurfaceCarrier::Plane {
            origin: points[0],
            normal: Vector3::z(),
        };
        let cylinder = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let result = path(1, false)
            .segment_analytic_trims(&points, &[plane, cylinder], options())
            .ok_or("segmentation")?;
        assert_eq!(result.pieces.len(), 1);
        let SharedTrimGeometry::Circle(circle) = result.pieces[0].geometry else {
            return Err("missing axial intersection".into());
        };
        assert!((circle.radius - 1.0).abs() < 1.0e-10);
        assert!((circle.signed_sweep_radians.abs() - std::f64::consts::FRAC_PI_2).abs() < 1.0e-10);
        let coincident = path(1, false)
            .segment_analytic_trims(&points, &[cylinder, cylinder], options())
            .ok_or("segmentation")?;
        assert!(matches!(
            coincident.pieces[0].geometry,
            SharedTrimGeometry::Unresolved
        ));
        Ok(())
    }

    #[test]
    fn two_samples_resolve_meridian_and_arbitrary_sphere_plane_sections()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(0.0, 0.0, 0.0);
        let torus = SurfaceCarrier::CircularRevolution {
            axial_center: center,
            axis: Vector3::z(),
            radial_center: 4.0,
            profile_radius: 1.0,
        };
        let points = [Point3::new(5.0, 0.0, 0.0), Point3::new(4.0, 0.0, 1.0)];
        let plane = SurfaceCarrier::Plane {
            origin: center,
            normal: Vector3::y(),
        };
        let result = path(1, false)
            .segment_analytic_trims(&points, &[torus, plane], options())
            .ok_or("segmentation")?;
        let SharedTrimGeometry::Circle(circle) = result.pieces[0].geometry else {
            return Err("missing meridian intersection".into());
        };
        assert!((circle.center - Point3::new(4.0, 0.0, 0.0)).norm() < 1.0e-10);
        let sphere = SurfaceCarrier::CircularRevolution {
            axial_center: center,
            axis: Vector3::z(),
            radial_center: 0.0,
            profile_radius: 1.0,
        };
        let radius = 0.75_f64.sqrt();
        let points = [Point3::new(0.5, radius, 0.0), Point3::new(0.5, 0.0, radius)];
        let plane = SurfaceCarrier::Plane {
            origin: points[0],
            normal: Vector3::x(),
        };
        let result = path(1, false)
            .segment_analytic_trims(&points, &[sphere, plane], options())
            .ok_or("segmentation")?;
        let SharedTrimGeometry::Circle(circle) = result.pieces[0].geometry else {
            return Err("missing sphere-plane intersection".into());
        };
        assert!((circle.center - Point3::new(0.5, 0.0, 0.0)).norm() < 1.0e-10);
        assert!((circle.radius - radius).abs() < 1.0e-10);
        Ok(())
    }

    #[test]
    fn spindle_meridian_trim_certifies_only_the_observed_radial_sheet()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(0.0, 0.0, 0.0);
        let torus = SurfaceCarrier::CircularRevolution {
            axial_center: center,
            axis: Vector3::z(),
            radial_center: 0.1,
            profile_radius: 1.0,
        };
        let cylinder = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.1, 0.0, 0.0),
            axis: Vector3::y(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let point = |angle: f64| {
            let (sin, cos) = angle.sin_cos();
            Point3::new(0.1 + cos, 0.0, sin)
        };
        // The circle has a negative-radial portion outside this observed
        // arc. The full supporting-circle proof must still reject it.
        let points = [point(0.0), point(0.4), point(0.8), point(1.2), point(1.5)];
        for ids in [vec![0, 1, 2, 3, 4], vec![4, 3, 2, 1, 0], vec![0, 1]] {
            let Some(SharedTrimGeometry::Circle(circle)) =
                accepted(&points, &ids, &[torus, cylinder], options())
            else {
                return Err("missing spindle meridian arc".into());
            };
            let fit = measure_path_curves(&points, &ids).ok_or("curve fit")?;
            let full = torus
                .measure(&crate::mesh_recovery::ContactCurveEvidence {
                    circle: Some(circle),
                    ..fit
                })
                .circle
                .ok_or("full circle")?;
            assert!(full.upper_bound > 0.1);
            assert!(
                torus
                    .measure_circle_arc(&circle, points[ids[0]])
                    .ok_or("finite arc")?
                    .upper_bound
                    < 1.0e-8
            );
        }
        // The same placed circle continued onto the opposite radial
        // sheet is not contained by the selected unsigned meridian.
        let crossing = [point(1.3), point(1.6), point(1.9), point(2.2)];
        let circle = super::super::curve_evidence::measure_placed_circle(
            &crossing,
            &[0, 1, 2, 3],
            Point3::new(0.1, 0.0, 0.0),
            Vector3::y(),
            1.0,
            1.0,
        )
        .ok_or("crossing circle")?;
        assert!(
            torus
                .measure_circle_arc(&circle, crossing[0])
                .ok_or("crossing arc")?
                .witnessed_max
                > 0.1
        );
        assert!(accepted(&crossing, &[0, 1, 2, 3], &[torus, cylinder], options()).is_none());
        Ok(())
    }

    #[test]
    fn open_path_three_samples_test_a_surface_supplied_circle()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 3.0_f64.sqrt() * 0.5, 0.0),
            Point3::new(-0.5, 3.0_f64.sqrt() * 0.5, 0.0),
        ];
        let cylinder = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let mut open = path(2, false);
        open.regions.1 = None;
        let fitted = open
            .segment_analytic_trims(&points, &[cylinder], options())
            .ok_or("trim")?;
        assert_eq!(fitted.pieces.len(), 1);
        assert!(matches!(
            fitted.pieces[0].geometry,
            SharedTrimGeometry::Circle(_)
        ));
        open.edges.truncate(1);
        let sparse = open
            .segment_analytic_trims(&points, &[cylinder], options())
            .ok_or("trim")?;
        assert!(matches!(
            sparse.pieces[0].geometry,
            SharedTrimGeometry::Unresolved
        ));
        Ok(())
    }

    #[test]
    fn sparse_three_point_circle_requires_independent_carrier_constraints()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
        ];
        let plane = SurfaceCarrier::Plane {
            origin: points[0],
            normal: Vector3::z(),
        };
        let cylinder = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let result = path(2, false)
            .segment_analytic_trims(&points, &[plane, cylinder], options())
            .ok_or("segmentation")?;
        assert_eq!(result.pieces.len(), 1);
        assert_eq!(result.pieces[0].edge_count, 2);
        assert!(matches!(
            result.pieces[0].geometry,
            SharedTrimGeometry::Circle(_)
        ));
        let planar = path(2, false)
            .segment_analytic_trims(&points, &[plane, plane], options())
            .ok_or("planar segmentation")?;
        assert!(
            planar
                .pieces
                .iter()
                .all(|piece| !matches!(piece.geometry, SharedTrimGeometry::Circle(_)))
        );
        Ok(())
    }
}
