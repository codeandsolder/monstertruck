use monstertruck_fillet::{FilletIntersectionCurve, ParameterCurveLinear};
use monstertruck_geometry::prelude::*;

use crate::{Curve, Curve2D, Surface};

impl TryFrom<Surface> for NurbsSurface<Vector4> {
    type Error = ();
    fn try_from(surface: Surface) -> std::result::Result<Self, ()> {
        match surface {
            Surface::Plane(plane) => Ok(NurbsSurface::from(BsplineSurface::from(plane))),
            Surface::BsplineSurface(bsp) => Ok(NurbsSurface::from(bsp)),
            Surface::NurbsSurface(ns) => Ok(ns),
            // The exact rational net -- the SAME one these faces WERE before
            // spec 012 U1.2 moved them onto the analytic variants, so the
            // fillet path sees no change.
            Surface::Sphere(_) | Surface::Torus(_) => surface
                .try_into_homogeneous_bspline_surface()
                .map(NurbsSurface::new)
                .ok_or(()),
            Surface::Revolution(_) | Surface::Tspline(_) => Err(()),
        }
    }
}
// From<NurbsSurface<Vector4>> for Surface -- provided by derive_more::From

impl TryFrom<Curve> for NurbsCurve<Vector4> {
    type Error = ();
    // Exact or refuse: boolean seam edges are `IntersectionCurve`s whose
    // exact leader must survive conversion; a sampled polyline stand-in would
    // silently degrade downstream fillet/chamfer geometry. Parameter curves
    // and non-exact leader chains have no exact NURBS form, so they refuse
    // here and surface as `FilletError::UnsupportedGeometry`.
    fn try_from(curve: Curve) -> std::result::Result<Self, ()> {
        curve
            .try_into_homogeneous_bspline_curve()
            .map(NurbsCurve::new)
            .ok_or(())
    }
}
// From<NurbsCurve<Vector4>> for Curve -- provided by derive_more::From

impl From<ParameterCurveLinear> for Curve {
    fn from(c: ParameterCurveLinear) -> Self {
        let (line, surface) = c.decompose();
        Curve::ParameterCurve(ParameterCurve::new(
            Curve2D::Line(line),
            Box::new(Surface::NurbsSurface(surface)),
        ))
    }
}

impl From<FilletIntersectionCurve> for Curve {
    fn from(curve: FilletIntersectionCurve) -> Self {
        let (surface0, surface1, leader) = curve.destruct();
        Curve::Intersection(SurfaceCurve::with_boundaries(
            Box::new(Surface::NurbsSurface(*surface0)),
            Box::new(Surface::NurbsSurface(*surface1)),
            Box::new(leader.into()),
            None,
            None,
        ))
    }
}

#[cfg(test)]
mod tests;
