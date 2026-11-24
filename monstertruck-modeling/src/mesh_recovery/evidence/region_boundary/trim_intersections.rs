//! Analytic circle proposals supplied by incident carrier geometry rather
//! than a three-point path fit. All proposals still require full containment
//! in both independently recovered incident surfaces.

use super::{ContactCircleEvidence, SurfaceCarrier, curve_evidence::measure_placed_circle};
use crate::{InnerSpace, Point3, Vector3};

#[derive(Clone, Copy)]
struct PlacedCircle {
    center: Point3,
    normal: Vector3,
    radius: f64,
}

fn surface_circles(
    carrier: SurfaceCarrier,
    other: SurfaceCarrier,
    sample: Point3,
    tolerance: f64,
) -> Vec<PlacedCircle> {
    match carrier {
        SurfaceCarrier::Plane { .. } => Vec::new(),
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } => {
            let axial = (sample - origin).dot(axis);
            vec![PlacedCircle {
                center: origin + axis * axial,
                normal: axis,
                radius: slope.mul_add(axial, radius_at_origin),
            }]
        }
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } => {
            let mut circles = Vec::new();
            let delta = sample - axial_center;
            let axial = delta.dot(axis);
            let offset = delta - axis * axial;
            // Constant-axial cross sections can have two distinct positive
            // radii. Both are proposed; the observed sheet and the other
            // surface decide which, if either, is the common curve.
            let squared = profile_radius.mul_add(profile_radius, -axial * axial);
            if squared >= 0.0 {
                let radial = squared.sqrt();
                for radius in [radial_center - radial, radial_center + radial] {
                    circles.push(PlacedCircle {
                        center: axial_center + axis * axial,
                        normal: axis,
                        radius,
                    });
                }
            }
            // A constant azimuth section is a meridian circle. If its
            // continuation crosses the axis onto a wrong radial branch,
            // the unsquared full-carrier containment check rejects it.
            if offset.magnitude() > tolerance {
                let direction = offset.normalize();
                circles.push(PlacedCircle {
                    center: axial_center + direction * radial_center,
                    normal: axis.cross(direction).normalize(),
                    radius: profile_radius,
                });
            }
            // A sphere is the zero-offset circular meridian; the arbitrary
            // plane section is independent of the chosen rotational axis.
            if radial_center.abs() <= tolerance
                && let SurfaceCarrier::Plane { origin, normal } = other
            {
                let distance = (axial_center - origin).dot(normal);
                let squared = profile_radius.mul_add(profile_radius, -distance * distance);
                if squared > 0.0 {
                    circles.push(PlacedCircle {
                        center: axial_center - normal * distance,
                        normal,
                        radius: squared.sqrt(),
                    });
                }
            }
            circles
        }
    }
}

pub(super) fn supplied_circles(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    tolerance: f64,
) -> Vec<ContactCircleEvidence> {
    let [first, second] = carriers else {
        return Vec::new();
    };
    let Some(&sample) = ids.first().and_then(|&id| points.get(id)) else {
        return Vec::new();
    };
    surface_circles(*first, *second, sample, tolerance)
        .into_iter()
        .chain(surface_circles(*second, *first, sample, tolerance))
        .filter_map(|circle| {
            measure_placed_circle(
                points,
                ids,
                circle.center,
                circle.normal,
                circle.radius,
                0.0,
            )
        })
        .filter(|circle| {
            circle.max_plane_residual <= tolerance && circle.max_vertex_radial_residual <= tolerance
        })
        .collect()
}

pub(super) fn same_circle(
    left: &ContactCircleEvidence,
    right: &ContactCircleEvidence,
    tolerance: f64,
) -> bool {
    (left.center - right.center).magnitude() <= tolerance
        && (left.radius - right.radius).abs() <= tolerance
        && left.normal.cross(right.normal).magnitude() * left.radius.max(right.radius) <= tolerance
}
