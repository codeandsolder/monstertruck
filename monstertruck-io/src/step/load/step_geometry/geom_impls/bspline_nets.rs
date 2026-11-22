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
    fn homogeneous_net_covers_domain(&self) -> bool { self.0.homogeneous_net_covers_domain() }

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
    fn homogeneous_net_covers_domain(&self) -> bool {
        match self {
            Self::Line(curve) => curve.homogeneous_net_covers_domain(),
            Self::Conic(Conic3D::Ellipse(curve)) => curve.homogeneous_net_covers_domain(),
            Self::BsplineCurve(curve) => curve.homogeneous_net_covers_domain(),
            Self::NurbsCurve(curve) => curve.homogeneous_net_covers_domain(),
            Self::SurfaceCurve(curve) => curve.leader().homogeneous_net_covers_domain(),
            Self::Conic(_)
            | Self::ParameterCurve(_)
            | Self::Polyline(_)
            | Self::IntersectionCurve(_) => false,
        }
    }

    fn try_into_homogeneous_bspline_curve(&self) -> Option<BsplineCurve<Vector4>> {
        match self {
            Self::Line(line) => line.try_into_homogeneous_bspline_curve(),
            Self::Conic(Conic3D::Ellipse(curve)) => curve.try_into_homogeneous_bspline_curve(),
            Self::BsplineCurve(curve) => curve.try_into_homogeneous_bspline_curve(),
            // A STEP surface curve evaluates its own leader directly.
            Self::SurfaceCurve(curve) => curve.leader().try_into_homogeneous_bspline_curve(),
            Self::NurbsCurve(curve) => curve.try_into_homogeneous_bspline_curve(),
            // A projected surface intersection evaluates Newton-corrected points;
            // its marching leader is not generally the same exact curve.
            Self::Conic(_)
            | Self::ParameterCurve(_)
            | Self::Polyline(_)
            | Self::IntersectionCurve(_) => None,
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

// Both spline conversions and the enclosure predicate must dispatch across
// the same set of STEP surface carriers, including their wrapped charts.
macro_rules! dispatch_surface_net {
    ($self:expr, $method:ident $(, $argument:expr)*) => {
        match $self {
            Surface::Elementary(ElementarySurface::Plane(surface)) => surface.$method($($argument),*),
            Surface::Elementary(ElementarySurface::Sphere(surface)) => surface.$method($($argument),*),
            Surface::Elementary(ElementarySurface::CylindricalSurface(surface) | ElementarySurface::ConicalSurface(surface)) => surface.$method($($argument),*),
            Surface::Elementary(ElementarySurface::ToroidalSurface(surface)) => surface.$method($($argument),*),
            Surface::Sweep(SweepSurface::ExtrusionSurface(surface)) => surface.$method($($argument),*),
            Surface::Sweep(SweepSurface::RevolutionSurface(surface)) => surface.$method($($argument),*),
            Surface::Bspline(surface) => surface.$method($($argument),*),
            Surface::Nurbs(surface) => surface.$method($($argument),*),
        }
    };
}

impl TryIntoBsplineSurface for Surface {
    fn try_into_bspline_surface(&self) -> Option<BsplineSurface<Point3>> {
        dispatch_surface_net!(self, try_into_bspline_surface)
    }
}

impl TryIntoHomogeneousBsplineSurface for Surface {
    fn homogeneous_net_covers_domain(&self) -> bool {
        dispatch_surface_net!(self, homogeneous_net_covers_domain)
    }

    fn try_into_homogeneous_bspline_surface(&self) -> Option<BsplineSurface<Vector4>> {
        dispatch_surface_net!(self, try_into_homogeneous_bspline_surface)
    }

    fn try_into_homogeneous_bspline_surface_over(
        &self,
        parameter_range: Option<SurfaceParameterRectangle>,
    ) -> Option<HomogeneousSurfaceConversion> {
        dispatch_surface_net!(
            self,
            try_into_homogeneous_bspline_surface_over,
            parameter_range
        )
    }
}

impl SupportsExactPatchDomains for Surface {
    fn supports_exact_patch_domains(&self) -> bool {
        matches!(self, Self::Bspline(_) | Self::Nurbs(_))
    }
}
