mod planar_network;
use planar_network::try_round_planar_network;

use monstertruck_fillet::{FilletIntersectionCurve, ParameterCurveLinear};
use monstertruck_geometry::prelude::*;

use crate::{Curve, Curve2D, Surface};

/// Fillet the native modeling shell, preferring simultaneous analytic
/// rolling-ball construction for compatible planar corner networks.
///
/// # Errors
/// Returns the generic fillet solver's error for other geometry or an
/// invalid request. No partial analytic network is ever committed.
pub fn fillet_edges(
    shell: &mut crate::Shell,
    edges: &[crate::Edge],
    options: Option<&monstertruck_fillet::FilletOptions>,
) -> std::result::Result<(), monstertruck_fillet::FilletError> {
    use monstertruck_fillet::{FilletOptions, FilletProfile, RadiusSpec};
    let default = FilletOptions::default();
    let options = options.unwrap_or(&default);
    if matches!(options.profile, FilletProfile::Round) {
        let constant;
        let radii = match &options.radius {
            RadiusSpec::Constant(radius) => {
                constant = vec![*radius; edges.len()];
                Some(constant.as_slice())
            }
            RadiusSpec::PerEdge(radii) => Some(radii.as_slice()),
            RadiusSpec::Variable(_) => None,
        };
        if let Some(rebuilt) = radii.and_then(|radii| try_round_planar_network(shell, edges, radii))
        {
            *shell = rebuilt;
            return Ok(());
        }
    }
    monstertruck_fillet::fillet_edges_generic(shell, edges, Some(options))
}

impl TryFrom<Surface> for NurbsSurface<Vector4> {
    type Error = ();
    fn try_from(surface: Surface) -> std::result::Result<Self, ()> {
        match surface {
            Surface::Plane(plane) => Ok(Self::from(BsplineSurface::from(plane))),
            Surface::BsplineSurface(bsp) => Ok(Self::from(bsp)),
            Surface::NurbsSurface(ns) => Ok(ns),
            // The exact rational net -- the SAME one these faces WERE before
            // spec 012 U1.2 moved them onto the analytic variants, so the
            // fillet path sees no change.
            Surface::Sphere(_) | Surface::Torus(_) => surface
                .try_into_homogeneous_bspline_surface()
                .map(Self::new)
                .ok_or(()),
            Surface::Revolution(_) | Surface::Tspline(_) => Err(()),
        }
    }
}
// From<NurbsSurface<Vector4>> for Surface -- provided by derive_more::From

impl TryFrom<Curve> for NurbsCurve<Vector4> {
    type Error = ();
    // Exact or refuse: surface intersections can evaluate to Newton-projected
    // points that differ from the marching B-spline leader. Such curves have
    // no certified exact NURBS conversion, and must be refused rather than
    // silently changing fillet/chamfer geometry. Parameter curves likewise
    // refuse unless an exact representation is independently available.
    fn try_from(curve: Curve) -> std::result::Result<Self, ()> {
        curve
            .try_into_homogeneous_bspline_curve()
            .map(Self::new)
            .ok_or(())
    }
}
// From<NurbsCurve<Vector4>> for Curve -- provided by derive_more::From

impl From<ParameterCurveLinear> for Curve {
    fn from(c: ParameterCurveLinear) -> Self {
        let (line, surface) = c.decompose();
        Self::ParameterCurve(ParameterCurve::new(
            Curve2D::Line(line),
            Box::new(Surface::NurbsSurface(surface)),
        ))
    }
}

impl From<FilletIntersectionCurve> for Curve {
    fn from(curve: FilletIntersectionCurve) -> Self {
        let (surface0, surface1, leader) = curve.destruct();
        Self::Intersection(SurfaceCurve::with_boundaries(
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
