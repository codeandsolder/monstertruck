use std::collections::HashMap;

use monstertruck_geometry::prelude::*;
use monstertruck_traits::{BoundedCurve, ParametricCurve, ParametricSurface};

use super::error::FilletError;
use super::types::{self, Curve, ParameterCurveLinear};

type InternalShell = types::Shell;
const SURFACE_SAMPLE_COUNT: usize = 64;

/// Intersection curve type used internally by fillet operations.
pub type FilletIntersectionCurve =
    IntersectionCurve<ParameterCurveLinear, Box<NurbsSurface<Vector4>>, Box<NurbsSurface<Vector4>>>;

/// Surface types that can participate in fillet operations.
///
/// Automatically implemented for any type satisfying the bounds.
pub trait FilletableSurface:
    Clone
    + ParametricSurface<Point = Point3>
    + TryInto<NurbsSurface<Vector4>>
    + From<NurbsSurface<Vector4>> {
    /// Converts this surface to a NURBS surface used by internal fillet logic.
    fn to_nurbs_surface(&self) -> Option<NurbsSurface<Vector4>> {
        self.clone()
            .try_into()
            .ok()
            .or_else(|| sample_surface_to_nurbs(self, SURFACE_SAMPLE_COUNT))
    }
}

impl<T> FilletableSurface for T where T: Clone
        + ParametricSurface<Point = Point3>
        + TryInto<NurbsSurface<Vector4>>
        + From<NurbsSurface<Vector4>>
{
}

/// Curve types that can participate in fillet operations.
///
/// Automatically implemented for any type satisfying the bounds.
pub trait FilletableCurve:
    Clone
    + ParametricCurve<Point = Point3>
    + BoundedCurve
    + TryInto<NurbsCurve<Vector4>>
    + From<NurbsCurve<Vector4>>
    + From<ParameterCurveLinear>
    + From<FilletIntersectionCurve> {
    /// Converts this curve to the exact NURBS curve used by internal fillet
    /// logic, or `None` when the curve has no exact NURBS representation.
    ///
    /// Refusal surfaces as [`FilletError::UnsupportedGeometry`] from
    /// [`fillet_edges_generic`](super::fillet_edges_generic); it must never
    /// be papered over with a sampled approximation, or exact boolean seam
    /// leaders would silently degrade into polylines.
    fn to_nurbs_curve(&self) -> Option<NurbsCurve<Vector4>> { self.clone().try_into().ok() }
}

impl<T> FilletableCurve for T where T: Clone
        + ParametricCurve<Point = Point3>
        + BoundedCurve
        + TryInto<NurbsCurve<Vector4>>
        + From<NurbsCurve<Vector4>>
        + From<ParameterCurveLinear>
        + From<FilletIntersectionCurve>
{
}

fn sample_surface_to_nurbs<S: ParametricSurface<Point = Point3>>(
    surface: &S,
    sample_count: usize,
) -> Option<NurbsSurface<Vector4>> {
    let (u_range, v_range) = surface.try_range_tuple();
    let ((u0, u1), (v0, v1)) = u_range.zip(v_range)?;
    let control_points: Vec<Vec<Point3>> = (0..=sample_count)
        .map(|iu| {
            let u =
                u0 + (u1 - u0) * (crate::count_as_f64(iu)) / (crate::count_as_f64(sample_count));
            (0..=sample_count)
                .map(|iv| {
                    let v = v0
                        + (v1 - v0) * (crate::count_as_f64(iv))
                            / (crate::count_as_f64(sample_count));
                    surface.evaluate(u, v)
                })
                .collect()
        })
        .collect();
    let u_knot = KnotVector::uniform_knot(1, sample_count);
    let v_knot = KnotVector::uniform_knot(1, sample_count);
    Some(NurbsSurface::from(BsplineSurface::new_unchecked(
        (u_knot, v_knot),
        control_points,
    )))
}

// TryFrom for extracting NurbsCurve from internal Curve type.
impl TryFrom<Curve> for NurbsCurve<Vector4> {
    type Error = ();
    fn try_from(curve: Curve) -> std::result::Result<Self, ()> {
        match curve {
            Curve::NurbsCurve(c) => Ok(c),
            _ => Err(()),
        }
    }
}

/// Convert an external shell to internal fillet types.
///
/// Returns the internal shell and the internal `EdgeId`s corresponding to
/// the selected external edges (matched by endpoint positions).
pub fn convert_shell_in<C: FilletableCurve, S: FilletableSurface>(
    shell: &monstertruck_topology::Shell<Point3, C, S>,
    edges: &[monstertruck_topology::Edge<Point3, C>],
) -> std::result::Result<(InternalShell, Vec<types::EdgeId>), FilletError> {
    let internal_shell: InternalShell = shell
        .try_mapped(
            |p| Some(*p),
            |c| c.to_nurbs_curve().map(Curve::NurbsCurve),
            FilletableSurface::to_nurbs_surface,
        )
        .ok_or(FilletError::UnsupportedGeometry {
            context: "failed to convert shell curves or surfaces to NURBS",
        })?;

    // Shell mapping retains face, boundary and edge order, even when several
    // distinct edges have identical endpoints. Match their topology identities,
    // not their coordinates (which are not unique in a CAD B-rep).
    let mut mapped_ids = HashMap::new();
    for (source_face, mapped_face) in shell.face_iter().zip(internal_shell.face_iter()) {
        for (source_wire, mapped_wire) in source_face
            .absolute_boundaries()
            .iter()
            .zip(mapped_face.absolute_boundaries())
        {
            for (source_edge, mapped_edge) in source_wire.iter().zip(mapped_wire.iter()) {
                mapped_ids.insert(source_edge.id(), mapped_edge.id());
            }
        }
    }
    let internal_edge_ids = edges
        .iter()
        .map(|edge| {
            mapped_ids
                .get(&edge.id())
                .copied()
                .ok_or(FilletError::EdgeNotFound)
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;

    Ok((internal_shell, internal_edge_ids))
}

/// Convert an internal fillet shell back to external types.
pub fn convert_shell_out<C: FilletableCurve, S: FilletableSurface>(
    shell: &InternalShell,
) -> std::result::Result<monstertruck_topology::Shell<Point3, C, S>, FilletError> {
    shell
        .try_mapped(
            |p| Some(*p),
            |c| {
                Some(match c {
                    Curve::NurbsCurve(nc) => C::from(nc.clone()),
                    Curve::ParameterCurve(pc) => C::from(pc.clone()),
                    Curve::IntersectionCurve(ic) => C::from(ic.clone()),
                })
            },
            |s| Some(S::from(s.clone())),
        )
        .ok_or(FilletError::UnsupportedGeometry {
            context: "failed to convert internal shell back to external types",
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_topology::{Edge, Face, Shell, Vertex, Wire};

    #[test]
    fn conversion_preserves_distinct_edges_with_identical_endpoints() {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let vertices = Vertex::from_points(points);
        let make_edge = |i: usize, j: usize| {
            let curve =
                BsplineCurve::new_unchecked(KnotVector::bezier_knot(1), vec![points[i], points[j]]);
            Edge::new_unchecked(
                &vertices[i],
                &vertices[j],
                Curve::NurbsCurve(NurbsCurve::from(curve)),
            )
        };
        let first = make_edge(0, 1);
        let second = make_edge(0, 1);
        let surface = || {
            let knots = KnotVector::bezier_knot(1);
            NurbsSurface::from(BsplineSurface::new_unchecked(
                (knots.clone(), knots),
                vec![vec![points[0], points[3]], vec![points[1], points[2]]],
            ))
        };
        let make_face = |edge: Edge<_, _>| {
            let wire: Wire<_, _> =
                vec![edge, make_edge(1, 2), make_edge(2, 3), make_edge(3, 0)].into();
            Face::new_unchecked(vec![wire], surface())
        };
        let shell: Shell<_, _, _> =
            vec![make_face(first.clone()), make_face(second.clone())].into();
        let (converted, ids) = convert_shell_in(&shell, &[second, first]).unwrap();
        assert_eq!(ids.len(), 2);
        assert_ne!(ids[0], ids[1]);
        assert!(converted.edge_iter().any(|edge| edge.id() == ids[0]));
        assert!(converted.edge_iter().any(|edge| edge.id() == ids[1]));
    }
}
