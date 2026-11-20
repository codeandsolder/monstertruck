//! Plane-separated material unions.  There can be no volume intersection
//! between the operands: only their oppositely facing planar footprints need
//! overlaying. The other faces are retained without SSI or volume classification.
use super::{
    AltCurve, AltRegion, Alternative, BooleanTolerance, Face, ShapeOpsCurve, ShapeOpsSurface,
    Shell, Solid,
    planar_union::{overlay_one_source, overlay_planar_subtract},
    weld_and_convert,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use monstertruck_geometry::prelude::*;

#[derive(Clone, Copy)]
struct Separator {
    origin: Point3,
    normal: Vector3,
}

impl Separator {
    fn from_face<C, S: ShapeOpsSurface>(face: &Face<Point3, C, S>) -> Option<Self> {
        let surface = face.oriented_surface();
        let AnalyticSurfaceKind::Plane(plane) = surface.try_into_analytic_surface_kind()? else {
            return None;
        };
        // A spline's carrier recognition is unoriented; its control-net
        // ordering need not agree with the regular chart normal.
        let point = face.vertex_iter().next()?.point();
        let (u, v) = surface.search_parameter(point, None, 100)?;
        let normal = surface.normal(u, v);
        let alignment = normal.dot(plane.normal());
        if !alignment.is_finite() || alignment.abs() < 1.0 - DEFAULT_NUMERICAL_TOLERANCES.angular {
            return None;
        }
        Some(Self {
            origin: plane.origin(),
            normal: plane.normal() * alignment.signum(),
        })
    }

    fn height(self, point: Point3) -> f64 { (point - self.origin).dot(self.normal) }
}

struct FaceBounds {
    points: Vec<Point3>,
    outward_plane: Option<Separator>,
}

/// Planes and curve extrusions have no extremum of a linear functional
/// away from their bounded wires: each point lies on a straight ruling
/// reaching the trim boundary. Positive rational weights certify that the
/// complete boundary curve lies inside its Euclidean control hull.
fn face_bounds<C: ShapeOpsCurve<S> + TryIntoHomogeneousBsplineCurve, S: ShapeOpsSurface>(
    face: &Face<Point3, C, S>,
) -> Option<FaceBounds> {
    if !face.surface().has_boundary_linear_extrema() {
        return None;
    }
    let mut points = Vec::new();
    for edge in face.edge_iter() {
        let curve = edge.curve().try_into_homogeneous_bspline_curve()?;
        for control in curve.control_points() {
            if !control.w.is_finite() || control.w <= 0.0 {
                return None;
            }
            let p = Point3::from_vec(control.truncate() / control.w);
            if ![p.x, p.y, p.z].into_iter().all(f64::is_finite) {
                return None;
            }
            points.push(p);
        }
    }
    (!points.is_empty()).then(|| FaceBounds {
        points,
        outward_plane: Separator::from_face(face),
    })
}

fn shell_bounds<C: ShapeOpsCurve<S> + TryIntoHomogeneousBsplineCurve, S: ShapeOpsSurface>(
    solid: &Solid<Point3, C, S>,
) -> Option<Vec<FaceBounds>> {
    solid.face_iter().map(face_bounds).collect()
}

/// Partition a closed boundary into caps on the separating plane and faces
/// strictly on the material side (possibly incident on the plane).
fn cap_partition(
    bounds: &[FaceBounds],
    separator: Separator,
    sign: f64,
    tol: f64,
) -> Option<Vec<bool>> {
    let mut caps = Vec::with_capacity(bounds.len());
    let mut has_interior = false;
    for face in bounds {
        let mut coplanar = true;
        for &point in &face.points {
            let distance = sign * separator.height(point);
            if !distance.is_finite() || distance > tol {
                return None;
            }
            coplanar &= distance.abs() <= tol;
            has_interior |= distance < -tol;
        }
        if coplanar {
            let plane = face.outward_plane?;
            // The direction is tested independently from the length budget.
            // An inward-facing or nonplanar contact cannot be discarded.
            if (plane.normal - separator.normal * sign).magnitude()
                > BooleanTolerance::new(tol).angular
                || separator.height(plane.origin).abs() > tol
            {
                return None;
            }
        }
        caps.push(coplanar);
    }
    (has_interior && caps.contains(&true)).then_some(caps)
}

struct PartitionedRegions<C, S> {
    retained: Vec<AltRegion<C, S>>,
    contact: Vec<AltRegion<C, S>>,
}

fn regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid: &Solid<Point3, C, S>,
    caps: &[bool],
    offset: usize,
) -> PartitionedRegions<C, S> {
    let mut retained = Vec::new();
    let mut contact = Vec::new();
    let mapped: Solid<Point3, AltCurve<C, S>, S> =
        solid.mapped(|p| *p, |c| Alternative::FirstType(c.clone()), Clone::clone);
    for (index, (face, &is_cap)) in mapped.face_iter().zip(caps).enumerate() {
        let region = AltRegion::new(
            face.absolute_boundaries().clone(),
            face.surface(),
            face.orientation(),
            offset + index,
            offset + index,
        );
        if is_cap {
            contact.push(region);
        } else {
            retained.push(region);
        }
    }
    PartitionedRegions { retained, contact }
}

/// Try to construct the union using a proven separating plane. A successful
/// result is assembled from the original exact side geometry and the two
/// exposed (positive minus opposite) portions of the contact footprints.
pub(in crate::transversal) fn union<
    C: ShapeOpsCurve<S> + TryIntoHomogeneousBsplineCurve,
    S: ShapeOpsSurface,
>(
    lhs: &Solid<Point3, C, S>,
    rhs: &Solid<Point3, C, S>,
    tol: f64,
) -> Option<Solid<Point3, C, S>> {
    let lhs_bounds = shell_bounds(lhs)?;
    let rhs_bounds = shell_bounds(rhs)?;
    for separator in lhs_bounds.iter().filter_map(|face| face.outward_plane) {
        let Some(left_caps) = cap_partition(&lhs_bounds, separator, 1.0, tol) else {
            continue;
        };
        let Some(right_caps) = cap_partition(&rhs_bounds, separator, -1.0, tol) else {
            continue;
        };
        let PartitionedRegions {
            retained: mut left,
            contact: left_contact,
        } = regions(lhs, &left_caps, 0);
        let PartitionedRegions {
            retained: right,
            contact: right_contact,
        } = regions(rhs, &right_caps, lhs_bounds.len());
        // The two differences are independently computed on their own
        // oriented charts, so both exposed cap normals remain authoritative.
        for (positive, cutter) in [
            (&left_contact, &right_contact),
            (&right_contact, &left_contact),
        ] {
            for cap in positive {
                left.extend(overlay_planar_subtract(
                    std::slice::from_ref(cap),
                    cutter,
                    tol,
                )?);
            }
        }
        left.extend(right);
        let result: Shell<Point3, C, S> =
            weld_and_convert(&left, "planar contact union", BooleanTolerance::new(tol)).ok()?;
        let components = result.connected_components();
        if components.is_empty() || components.iter().any(|s| s.check_solid_boundary().is_err()) {
            return None;
        }
        return Solid::new(components).ok();
    }
    None
}

/// Union co-oriented planar material regions, including overlapping profiles
/// and disjoint islands. Every boundary's complete curve hull must lie on
/// the common plane. Use a fresh plane chart rather than extending the
/// finite patch of any one input face.
pub(in crate::transversal) fn profiles<C, S>(
    faces: &[Face<Point3, C, S>],
    tol: f64,
) -> Option<Shell<Point3, C, S>>
where
    C: ShapeOpsCurve<S> + TryIntoHomogeneousBsplineCurve,
    S: ShapeOpsSurface,
    Plane: ToSameGeometry<S>,
{
    let reference = Separator::from_face(faces.first()?)?;
    let surface = faces.first()?.oriented_surface();
    let AnalyticSurfaceKind::Plane(mut plane) = surface.try_into_analytic_surface_kind()? else {
        return None;
    };
    if plane.normal().dot(reference.normal) < 0.0 {
        plane.invert();
    }
    let chart = plane.to_same_geometry();
    let mut regions = Vec::with_capacity(faces.len());
    for (i, face) in faces.iter().enumerate() {
        let bounds = face_bounds(face)?;
        let outward = bounds.outward_plane?;
        if (outward.normal - reference.normal).magnitude() > DEFAULT_NUMERICAL_TOLERANCES.angular
            || bounds
                .points
                .iter()
                .any(|&p| reference.height(p).abs() > tol)
        {
            return None;
        }
        let mapped: Face<Point3, AltCurve<C, S>, S> = face.mapped(
            |p| *p,
            |curve| Alternative::FirstType(curve.clone()),
            Clone::clone,
        );
        regions.push(AltRegion::new(
            mapped.boundaries(),
            chart.clone(),
            true,
            i,
            i,
        ));
    }
    let union = overlay_one_source(&regions, tol)?;
    if union.is_empty() {
        return None;
    }
    weld_and_convert(&union, "planar profile union", BooleanTolerance::new(tol)).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_modeling::{Curve, Surface, builder};
    type ModelSolid = Solid<Point3, Curve, Surface>;

    fn box_at(origin: Point3, size: Vector3) -> ModelSolid {
        let v = builder::vertex(origin);
        let e = builder::extrude(&v, Vector3::unit_x() * size.x);
        let f = builder::extrude(&e, Vector3::unit_y() * size.y);
        builder::extrude(&f, Vector3::unit_z() * size.z)
    }

    #[test]
    fn planar_contact_union_nodes_crossing_caps_and_refuses_volume_overlap() {
        let bottom = box_at(Point3::origin(), Vector3::new(5., 3., 2.));
        let top = box_at(Point3::new(2., -1., 2.), Vector3::new(2., 5., 1.));
        let result = union(&bottom, &top, 1e-6).expect("crossing cap footprints");
        assert_eq!(result.boundaries().len(), 1);
        assert!(result.is_geometric_consistent());
        assert!(result.boundaries()[0].check_solid_boundary().is_ok());
        let reverse = union(&top, &bottom, 1e-6).expect("reversed operands");
        assert_eq!(reverse.boundaries().len(), 1);
        assert!(reverse.is_geometric_consistent());
        let shifted = builder::translated(&top, Vector3::new(0., 0., -0.1));
        assert!(union(&bottom, &shifted, 1e-6).is_none());
        let gap = builder::translated(&top, Vector3::new(0., 0., 0.1));
        assert!(union(&bottom, &gap, 1e-6).is_none());
    }

    #[test]
    fn planar_profile_union_preserves_nested_circular_islands() {
        use monstertruck_modeling::{primitive, profile};
        let circle = |r: f64| {
            primitive::circle::<Curve>(
                Point3::new(r, 0., 0.),
                Point3::origin(),
                Vector3::unit_z(),
                4,
            )
        };
        let annulus: Face<Point3, Curve, Surface> =
            profile::attach_plane_normalized(vec![circle(3.), circle(2.)]).unwrap();
        let disk: Face<Point3, Curve, Surface> =
            profile::attach_plane_normalized(vec![circle(1.)]).unwrap();
        let regions = profiles(&[annulus, disk], 1e-6).expect("annulus and inner island");
        assert_eq!(regions.len(), 2);
        let mut counts: Vec<_> = regions.iter().map(|face| face.boundaries().len()).collect();
        counts.sort_unstable();
        assert_eq!(counts, [1, 2]);
        let mut boundaries = Vec::new();
        for face in regions.iter() {
            let component: ModelSolid = builder::extrude(face, Vector3::unit_z());
            boundaries.extend(component.into_boundaries());
        }
        let lower = Solid::new(boundaries).unwrap();
        assert_eq!(lower.boundaries().len(), 2);
        let cap: Face<Point3, Curve, Surface> =
            profile::attach_plane_normalized(vec![circle(3.)]).unwrap();
        let upper: ModelSolid = builder::extrude(&cap, Vector3::unit_z());
        let upper = builder::translated(&upper, Vector3::unit_z());
        let joined = union(&lower, &upper, 1e-6).expect("upper slab connects both islands");
        assert_eq!(joined.boundaries().len(), 1);
        assert!(joined.is_geometric_consistent());
    }
}
