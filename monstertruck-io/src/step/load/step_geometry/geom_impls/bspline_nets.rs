//! Flattening STEP geometry onto B-spline and rational (homogeneous)
//! B-spline nets, and the predicate for which surfaces carry an exact
//! patch domain.

use super::{
    BsplineCurve, BsplineSurface, Conic3D, Curve3D, ElementarySurface,
    HomogeneousSurfaceConversion, Point3, Sphere, SupportsExactPatchDomains, Surface,
    SurfaceParameterRectangle, SweepSurface, TryIntoBsplineSurface, TryIntoHomogeneousBsplineCurve,
    TryIntoHomogeneousBsplineSurface, Vector4,
};

impl TryIntoHomogeneousBsplineSurface for Sphere {
    fn try_into_homogeneous_bspline_surface(&self) -> Option<BsplineSurface<Vector4>> {
        self.0.try_into_homogeneous_bspline_surface()
    }

    fn try_into_homogeneous_bspline_surface_over(
        &self,
        parameter_range: Option<SurfaceParameterRectangle>,
    ) -> Option<HomogeneousSurfaceConversion> {
        self.0
            .try_into_homogeneous_bspline_surface_over(parameter_range)
    }
}

impl TryIntoBsplineSurface for Sphere {
    fn try_into_bspline_surface(&self) -> Option<BsplineSurface<Point3>> {
        self.0.try_into_bspline_surface()
    }
}

impl TryIntoHomogeneousBsplineCurve for Curve3D {
    fn try_into_homogeneous_bspline_curve(&self) -> Option<BsplineCurve<Vector4>> {
        match self {
            Self::Line(line) => line.try_into_homogeneous_bspline_curve(),
            Self::Conic(Conic3D::Ellipse(curve)) => curve.try_into_homogeneous_bspline_curve(),
            Self::Conic(_) | Self::ParameterCurve(_) | Self::Polyline(_) => None,
            Self::BsplineCurve(curve) => curve.try_into_homogeneous_bspline_curve(),
            Self::SurfaceCurve(curve) => curve.leader().try_into_homogeneous_bspline_curve(),
            Self::IntersectionCurve(curve) => curve.leader().try_into_homogeneous_bspline_curve(),
            Self::NurbsCurve(curve) => curve.try_into_homogeneous_bspline_curve(),
        }
    }

    fn try_into_homogeneous_bspline_curve_over(
        &self,
        range: (f64, f64),
    ) -> Option<BsplineCurve<Vector4>> {
        match self {
            // Only a line has an exact analytic continuation past its own range;
            // every other variant keeps the trait's refusing default.
            Self::Line(line) => line.try_into_homogeneous_bspline_curve_over(range),
            _ => None,
        }
    }
}

impl TryIntoBsplineSurface for Surface {
    fn try_into_bspline_surface(&self) -> Option<BsplineSurface<Point3>> {
        match self {
            Self::Elementary(ElementarySurface::Plane(surface)) => {
                surface.try_into_bspline_surface()
            }
            Self::Elementary(ElementarySurface::Sphere(surface)) => {
                surface.try_into_bspline_surface()
            }
            Self::Elementary(
                ElementarySurface::CylindricalSurface(surface)
                | ElementarySurface::ConicalSurface(surface),
            ) => surface.try_into_bspline_surface(),
            Self::Elementary(ElementarySurface::ToroidalSurface(surface)) => {
                surface.try_into_bspline_surface()
            }
            Self::Sweep(SweepSurface::ExtrusionSurface(surface)) => {
                surface.try_into_bspline_surface()
            }
            Self::Sweep(SweepSurface::RevolutionSurface(surface)) => {
                surface.try_into_bspline_surface()
            }
            Self::Bspline(surface) => surface.try_into_bspline_surface(),
            Self::Nurbs(surface) => surface.try_into_bspline_surface(),
        }
    }
}

impl TryIntoHomogeneousBsplineSurface for Surface {
    fn try_into_homogeneous_bspline_surface(&self) -> Option<BsplineSurface<Vector4>> {
        match self {
            Self::Elementary(ElementarySurface::Plane(surface)) => {
                surface.try_into_homogeneous_bspline_surface()
            }
            Self::Elementary(ElementarySurface::Sphere(surface)) => {
                surface.try_into_homogeneous_bspline_surface()
            }
            Self::Elementary(
                ElementarySurface::CylindricalSurface(surface)
                | ElementarySurface::ConicalSurface(surface),
            ) => surface.try_into_homogeneous_bspline_surface(),
            Self::Elementary(ElementarySurface::ToroidalSurface(surface)) => {
                surface.try_into_homogeneous_bspline_surface()
            }
            Self::Sweep(SweepSurface::ExtrusionSurface(surface)) => {
                surface.try_into_homogeneous_bspline_surface()
            }
            Self::Sweep(SweepSurface::RevolutionSurface(surface)) => {
                surface.try_into_homogeneous_bspline_surface()
            }
            Self::Bspline(surface) => surface.try_into_homogeneous_bspline_surface(),
            Self::Nurbs(surface) => surface.try_into_homogeneous_bspline_surface(),
        }
    }

    fn try_into_homogeneous_bspline_surface_over(
        &self,
        parameter_range: Option<SurfaceParameterRectangle>,
    ) -> Option<HomogeneousSurfaceConversion> {
        match self {
            Self::Elementary(ElementarySurface::Plane(surface)) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
            Self::Elementary(ElementarySurface::Sphere(surface)) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
            Self::Elementary(
                ElementarySurface::CylindricalSurface(surface)
                | ElementarySurface::ConicalSurface(surface),
            ) => surface.try_into_homogeneous_bspline_surface_over(parameter_range),
            Self::Elementary(ElementarySurface::ToroidalSurface(surface)) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
            Self::Sweep(SweepSurface::ExtrusionSurface(surface)) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
            Self::Sweep(SweepSurface::RevolutionSurface(surface)) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
            Self::Bspline(surface) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
            Self::Nurbs(surface) => {
                surface.try_into_homogeneous_bspline_surface_over(parameter_range)
            }
        }
    }
}

impl SupportsExactPatchDomains for Surface {
    fn supports_exact_patch_domains(&self) -> bool {
        matches!(self, Self::Bspline(_) | Self::Nurbs(_))
    }
}
