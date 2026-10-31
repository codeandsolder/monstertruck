//! Curve behaviour on the STEP curve enums: evaluation, division, cutting,
//! inversion, endpoint snapping, parameter search and transformation.

use super::{
    BoundedCurve, BsplineCurve, Curve3D, CurveParameter, Cut, Instant, IntersectionCurve,
    Invertible, Matrix4, ParameterDivision1D, ParameterRange, ParametricCurve, Point3,
    SearchNearestParameter, SearchParameter, SearchParameterHint1D, SnapCurveEndpoints, Surface,
    SurfaceCurve3D, SurfaceCurveAssociatedGeometry, ToSameGeometry, Transformed, Vector3, env,
};

impl ParametricCurve for SurfaceCurve3D {
    type Point = Point3;
    type Vector = Vector3;

    fn evaluate(&self, t: f64) -> Self::Point { self.leader().evaluate(t) }

    fn derivative(&self, t: f64) -> Self::Vector { self.leader().derivative(t) }

    fn derivative_2(&self, t: f64) -> Self::Vector { self.leader().derivative_2(t) }

    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector { self.leader().derivative_n(n, t) }

    fn parameter_range(&self) -> ParameterRange { self.leader().parameter_range() }

    fn period(&self) -> Option<f64> { self.leader().period() }
}

impl BoundedCurve for SurfaceCurve3D {
    fn range_tuple(&self) -> (f64, f64) { self.leader().range_tuple() }
}

impl ParameterDivision1D for SurfaceCurve3D {
    type Point = Point3;

    fn try_parameter_division(
        &self,
        range: (f64, f64),
        tol: f64,
    ) -> Option<(Vec<f64>, Vec<Self::Point>)> {
        self.leader().try_parameter_division(range, tol)
    }

    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        self.leader().parameter_division(range, tol)
    }
}

impl Cut for SurfaceCurve3D {
    fn cut(&mut self, t: f64) -> Self {
        let leader = Box::new(self.leader_mut().cut(t));
        let associated_geometry = self
            .associated_geometry
            .iter_mut()
            .map(|entry| entry.split_at(t))
            .collect();
        Self::new(
            self.kind(),
            leader,
            associated_geometry,
            self.master_representation(),
        )
    }
}

impl SnapCurveEndpoints for SurfaceCurve3D {
    fn snap_endpoints(&mut self, front: Point3, back: Point3) {
        self.leader_mut().snap_endpoints(front, back);
    }
}

impl SnapCurveEndpoints for Curve3D {
    fn snap_endpoints(&mut self, front: Point3, back: Point3) {
        match self {
            Self::Polyline(curve) => curve.snap_endpoints(front, back),
            Self::SurfaceCurve(curve) => curve.snap_endpoints(front, back),
            Self::IntersectionCurve(curve) => curve.snap_endpoints(front, back),
            Self::Line(_)
            | Self::Conic(_)
            | Self::BsplineCurve(_)
            | Self::ParameterCurve(_)
            | Self::NurbsCurve(_) => {}
        }
    }
}

impl Invertible for SurfaceCurveAssociatedGeometry {
    fn invert(&mut self) {
        if let Self::ParameterCurve(curve) = self {
            curve.invert();
        }
    }
}

impl Invertible for SurfaceCurve3D {
    fn invert(&mut self) {
        self.leader_mut().invert();
        self.associated_geometry
            .iter_mut()
            .for_each(Invertible::invert);
    }
}

impl SearchParameter<CurveParameter> for SurfaceCurve3D {
    type Point = Point3;

    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        self.leader().search_parameter(point, hint, trials)
    }
}

impl SearchNearestParameter<CurveParameter> for SurfaceCurve3D {
    type Point = Point3;

    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        self.leader().search_nearest_parameter(point, hint, trials)
    }
}

impl Transformed<Matrix4> for SurfaceCurveAssociatedGeometry {
    fn transform_by(&mut self, trans: Matrix4) {
        match self {
            Self::ParameterCurve(curve) => curve.transform_by(trans),
            Self::Surface(surface) => surface.transform_by(trans),
        }
    }
}

impl Transformed<Matrix4> for SurfaceCurve3D {
    fn transform_by(&mut self, trans: Matrix4) {
        self.leader_mut().transform_by(trans);
        self.associated_geometry
            .iter_mut()
            .for_each(|entry| entry.transform_by(trans));
    }
}

impl ParameterDivision1D for Curve3D {
    type Point = Point3;

    fn try_parameter_division(
        &self,
        range: (f64, f64),
        tol: f64,
    ) -> Option<(Vec<f64>, Vec<Self::Point>)> {
        let debug_profile = env::var("MT_PROFILE_CURVE_DIVISION").is_ok();
        // Only consult the clock when actually profiling -- `Instant::now()`
        // panics on `wasm32-unknown-unknown` ("time not implemented"), so
        // an unconditional call here breaks browser STEP loading.
        let started = debug_profile.then(Instant::now);
        let result = match self {
            Self::Line(curve) => curve.try_parameter_division(range, tol),
            Self::Polyline(curve) => curve.try_parameter_division(range, tol),
            Self::Conic(curve) => curve.try_parameter_division(range, tol),
            Self::BsplineCurve(curve) => curve.try_parameter_division(range, tol),
            Self::ParameterCurve(curve) => curve.try_parameter_division(range, tol),
            Self::SurfaceCurve(curve) => curve.try_parameter_division(range, tol),
            Self::IntersectionCurve(curve) => curve.leader().try_parameter_division(range, tol),
            Self::NurbsCurve(curve) => curve.try_parameter_division(range, tol),
        };
        if let Some(started) = started {
            let kind = match self {
                Self::Line(_) => "Line",
                Self::Polyline(_) => "Polyline",
                Self::Conic(_) => "Conic",
                Self::BsplineCurve(_) => "BsplineCurve",
                Self::ParameterCurve(_) => "StepParameterCurve",
                Self::SurfaceCurve(_) => "SurfaceCurve",
                Self::IntersectionCurve(_) => "IntersectionCurve",
                Self::NurbsCurve(_) => "NurbsCurve",
            };
            eprintln!(
                "trace bool curve_division kind={kind} points={} tol={tol} elapsed_ms={:.3}",
                result.as_ref().map_or(0, |(_, points)| points.len()),
                started.elapsed().as_secs_f64() * 1000.0,
            );
        }
        result
    }

    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        let debug_profile = env::var("MT_PROFILE_CURVE_DIVISION").is_ok();
        // Same wasm-safety guard as `try_parameter_division` above.
        let started = debug_profile.then(Instant::now);
        let result = match self {
            Self::Line(curve) => curve.parameter_division(range, tol),
            Self::Polyline(curve) => curve.parameter_division(range, tol),
            Self::Conic(curve) => curve.parameter_division(range, tol),
            Self::BsplineCurve(curve) => curve.parameter_division(range, tol),
            Self::ParameterCurve(curve) => curve.parameter_division(range, tol),
            Self::SurfaceCurve(curve) => curve.parameter_division(range, tol),
            Self::IntersectionCurve(curve) => curve.leader().parameter_division(range, tol),
            Self::NurbsCurve(curve) => curve.parameter_division(range, tol),
        };
        if let Some(started) = started {
            let kind = match self {
                Self::Line(_) => "Line",
                Self::Polyline(_) => "Polyline",
                Self::Conic(_) => "Conic",
                Self::BsplineCurve(_) => "BsplineCurve",
                Self::ParameterCurve(_) => "StepParameterCurve",
                Self::SurfaceCurve(_) => "SurfaceCurve",
                Self::IntersectionCurve(_) => "IntersectionCurve",
                Self::NurbsCurve(_) => "NurbsCurve",
            };
            eprintln!(
                "trace bool curve_division kind={kind} points={} tol={tol} elapsed_ms={:.3}",
                result.1.len(),
                started.elapsed().as_secs_f64() * 1000.0,
            );
        }
        result
    }
}

impl ToSameGeometry<Curve3D> for SurfaceCurve3D {
    fn to_same_geometry(&self) -> Curve3D { Curve3D::SurfaceCurve(self.clone()) }
}

impl From<IntersectionCurve<BsplineCurve<Point3>, Surface, Surface>> for Curve3D {
    fn from(ic: IntersectionCurve<BsplineCurve<Point3>, Surface, Surface>) -> Self {
        let (surface0, surface1, leader) = ic.destruct();
        Self::IntersectionCurve(IntersectionCurve::new(
            Box::new(surface0),
            Box::new(surface1),
            Box::new(Self::BsplineCurve(leader)),
        ))
    }
}
