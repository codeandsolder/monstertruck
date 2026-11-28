use super::{
    ConstStepLength, Debug, ERR, FloatDisplay, Formatter, IndexSliceDisplay, Result, SliceDisplay,
    StepCurve, StepDisplay, StepFormat, StepLength, StepPcurve, StepSurface,
    canonical_step_encoding,
};
use monstertruck_geometry::prelude::*;
use monstertruck_mesh::PolylineCurve;
use monstertruck_modeling::{
    Conic2D as ModelingConic2D, Curve as ModelingCurve, Curve2D as ModelingCurve2D,
    Surface as ModelingSurface,
};

impl StepFormat for Point2 {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        f.write_fmt(format_args!(
            "#{idx} = CARTESIAN_POINT('', {coordinates});\n",
            coordinates = SliceDisplay(AsRef::<[f64; 2]>::as_ref(self)),
        ))
    }
}
impl_const_step_length!(Point2, 1);

impl StepFormat for Point3 {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        f.write_fmt(format_args!(
            "#{idx} = CARTESIAN_POINT('', {coordinates});\n",
            coordinates = SliceDisplay(AsRef::<[f64; 3]>::as_ref(self)),
        ))
    }
}
impl_const_step_length!(Point3, 1);

/// class for display `DIRECTION`.
#[derive(Clone, Debug, Copy)]
pub struct VectorAsDirection<V>(pub V);

impl StepFormat for VectorAsDirection<Vector2> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        f.write_fmt(format_args!(
            "#{idx} = DIRECTION('', {direction_ratios});\n",
            direction_ratios = SliceDisplay(AsRef::<[f64; 2]>::as_ref(&self.0)),
        ))
    }
}

impl StepFormat for VectorAsDirection<Vector3> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        f.write_fmt(format_args!(
            "#{idx} = DIRECTION('', {direction_ratios});\n",
            direction_ratios = SliceDisplay(AsRef::<[f64; 3]>::as_ref(&self.0)),
        ))
    }
}
impl_const_step_length!(VectorAsDirection<V>, 1, <V>);

impl<V> StepFormat for V
where
    V: InnerProduct<Scalar = f64>,
    VectorAsDirection<V>: StepFormat,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let magnitude = self.norm();
        let direction_idx = idx + 1;
        f.write_fmt(format_args!(
            "#{idx} = VECTOR('', #{direction_idx}, {magnitude});\n{direction}",
            direction = StepDisplay::new(VectorAsDirection(*self / magnitude), direction_idx),
            magnitude = FloatDisplay(magnitude),
        ))
    }
}
impl_const_step_length!(Vector2, 2);
impl_const_step_length!(Vector3, 2);

/// Wrapper that displays a transform matrix as a STEP `AXIS2_PLACEMENT_*` entity.
///
/// STEP encodes an oriented coordinate frame -- a placement -- as a
/// location point plus one (`AXIS2_PLACEMENT_2D`) or two
/// (`AXIS2_PLACEMENT_3D`) direction vectors. This wrapper takes a
/// native column-major transform matrix and emits the four
/// referenced entities: the placement itself, the location point, and
/// the direction(s).
#[derive(Clone, Copy, Debug)]
pub struct MatrixAsAxis<M>(pub M);

impl StepFormat for MatrixAsAxis<Matrix3> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let location_idx = idx + 1;
        let ref_direction_idx = idx + 2;
        let location = self.0.column(2).into_owned().to_point();
        let ref_direction = VectorAsDirection(self.0.column(0).into_owned().spatial());
        f.write_fmt(format_args!(
            "#{idx} = AXIS2_PLACEMENT_2D('', #{location_idx}, #{ref_direction_idx});\n",
        ))?;
        StepFormat::fmt(&location, location_idx, f)?;
        StepFormat::fmt(&ref_direction, ref_direction_idx, f)
    }
}
impl_const_step_length!(MatrixAsAxis<Matrix3>, 3);

impl StepFormat for MatrixAsAxis<Matrix4> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let location_idx = idx + 1;
        let axis_idx = idx + 2;
        let ref_direction_idx = idx + 3;
        let location = self.0.column(3).into_owned().to_point();
        let axis = VectorAsDirection(self.0.column(2).into_owned().spatial());
        let ref_direction = VectorAsDirection(self.0.column(0).into_owned().spatial());
        f.write_fmt(format_args!(
            "#{idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});\n",
        ))?;
        StepFormat::fmt(&location, location_idx, f)?;
        StepFormat::fmt(&axis, axis_idx, f)?;
        StepFormat::fmt(&ref_direction, ref_direction_idx, f)
    }
}
impl_const_step_length!(MatrixAsAxis<Matrix4>, 4);

impl<P> StepFormat for Line<P>
where
    P: AffinePoint + ConstStepLength + StepFormat,
    P::Diff: StepFormat,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let pnt_idx = idx + 1;
        let dir_idx = idx + 1 + P::LENGTH;
        f.write_fmt(format_args!(
            "#{idx} = LINE('', #{pnt_idx}, #{dir_idx});\n{pnt}{dir}",
            pnt = StepDisplay::new(self.0, pnt_idx),
            dir = StepDisplay::new(self.1 - self.0, dir_idx),
        ))
    }
}

impl<P> StepLength for Line<P>
where
    P: AffinePoint + ConstStepLength,
    P::Diff: ConstStepLength,
{
    fn step_length(&self) -> usize { <Self as ConstStepLength>::LENGTH }
}

impl<P> ConstStepLength for Line<P>
where
    P: AffinePoint + ConstStepLength,
    P::Diff: ConstStepLength,
{
    const LENGTH: usize = 1 + P::LENGTH + P::Diff::LENGTH;
}

impl<P> StepCurve for Line<P> {}

impl<P> StepFormat for PolylineCurve<P>
where P: Copy + ConstStepLength + StepFormat
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        f.write_fmt(format_args!(
            "#{idx} = POLYLINE('', {range});\n",
            range = IndexSliceDisplay(idx + 1..=idx + self.0.len())
        ))?;
        let closure = |(i, p): (usize, &P)| p.fmt(idx + 1 + i * P::LENGTH, f);
        self.0.iter().enumerate().try_for_each(closure)
    }
}

impl<P: ConstStepLength> StepLength for PolylineCurve<P> {
    fn step_length(&self) -> usize { 1 + self.0.len() * P::LENGTH }
}

impl<P> StepCurve for PolylineCurve<P> {}

impl<P> StepFormat for BsplineCurve<P>
where P: Copy + ConstStepLength + StepFormat
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let (knots, multi) = self.knot_vector().to_single_multi();
        let control_points_instances = self
            .control_points()
            .iter()
            .enumerate()
            .map(|(i, p)| StepDisplay::new(*p, idx + 1 + i * P::LENGTH))
            .collect::<Vec<_>>();
        f.write_fmt(format_args!(
            "#{idx} = B_SPLINE_CURVE_WITH_KNOTS('', {degree}, {control_points_list}, .UNSPECIFIED., .U., .U., {knot_multiplicities}, {knots}, .UNSPECIFIED.);\n{control_points_instances}",
            degree = self.degree(),
            control_points_list = IndexSliceDisplay((idx + 1..=idx + self.control_points().len() * P::LENGTH).step_by(P::LENGTH)),
			knot_multiplicities = SliceDisplay(&multi),
            knots = SliceDisplay(&knots),
            control_points_instances = SliceDisplay(&control_points_instances),
		))
    }
}

impl<P> StepLength for BsplineCurve<P> {
    fn step_length(&self) -> usize { self.control_points().len() + 1 }
}

impl<P> StepCurve for BsplineCurve<P> {}

impl<V> StepFormat for NurbsCurve<V>
where
    V: Homogeneous<Scalar = f64>,
    V::Point: ConstStepLength + StepFormat,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let (knots, multi) = self.knot_vector().to_single_multi();
        let control_points_instances = self
            .control_points()
            .iter()
            .enumerate()
            .map(|(i, v)| StepDisplay::new(v.to_point(), idx + 1 + i * V::Point::LENGTH))
            .collect::<Vec<_>>();
        let weights = self
            .control_points()
            .iter()
            .map(|v| v.weight())
            .collect::<Vec<_>>();
        f.write_fmt(format_args!(
            "#{idx} = (
    BOUNDED_CURVE()
    B_SPLINE_CURVE({degree}, {control_points_list}, .UNSPECIFIED., .U., .U.)
    B_SPLINE_CURVE_WITH_KNOTS({knot_multiplicities}, {knots}, .UNSPECIFIED.)
    CURVE()
    GEOMETRIC_REPRESENTATION_ITEM()
    RATIONAL_B_SPLINE_CURVE({weights})
    REPRESENTATION_ITEM('')
);\n{control_points_instances}",
            degree = self.degree(),
            control_points_list = IndexSliceDisplay(
                (idx + 1..=idx + self.control_points().len() * V::Point::LENGTH)
                    .step_by(V::Point::LENGTH)
            ),
            knot_multiplicities = SliceDisplay(&multi),
            knots = SliceDisplay(&knots),
            weights = SliceDisplay(&weights),
            control_points_instances = SliceDisplay(&control_points_instances),
        ))
    }
}

impl<V> StepLength for NurbsCurve<V> {
    fn step_length(&self) -> usize { self.control_points().len() + 1 }
}

impl<V> StepCurve for NurbsCurve<V> {}

impl StepFormat for Processor<TrimmedCurve<UnitCircle<Point2>>, Matrix3> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let transform = *self.transform();
        let position_idx = idx + 1;
        let location_idx = idx + 2;
        let ref_direction_idx = idx + 3;
        let r0 = transform.column(0).into_owned().norm();
        let r1 = transform.column(1).into_owned().norm();
        let ref_direction = VectorAsDirection(transform.column(0).into_owned().spatial() / r0);
        let location = transform.column(2).into_owned().to_point();
        if r0.near(&r1) {
            let r = FloatDisplay(r0);
            f.write_fmt(format_args!("#{idx} = CIRCLE('', #{position_idx}, {r});\n"))?;
        } else {
            let (r0, r1) = (FloatDisplay(r0), FloatDisplay(r1));
            f.write_fmt(format_args!(
                "#{idx} = ELLIPSE('', #{position_idx}, {r0}, {r1});\n"
            ))?;
        }
        f.write_fmt(format_args!(
            "#{position_idx} = AXIS2_PLACEMENT_2D('', #{location_idx}, #{ref_direction_idx});\n",
        ))?;
        StepFormat::fmt(&location, location_idx, f)?;
        StepFormat::fmt(&ref_direction, ref_direction_idx, f)
    }
}
impl_const_step_length!(Processor<TrimmedCurve<UnitCircle<Point2>>, Matrix3>, 4);

fn write_conic_placement_3d(transform: &Matrix4, idx: usize, f: &mut Formatter<'_>) -> Result {
    let (location_idx, axis_idx, ref_direction_idx) = (idx + 1, idx + 2, idx + 3);
    let location = transform.column(3).into_owned().to_point();
    let axis = VectorAsDirection(transform.column(2).into_owned().spatial().normalize());
    let direction =
        VectorAsDirection(transform.column(0).into_owned().spatial() / transform.column(0).norm());
    f.write_fmt(format_args!(
        "#{idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});\n"
    ))?;
    StepFormat::fmt(&location, location_idx, f)?;
    StepFormat::fmt(&axis, axis_idx, f)?;
    StepFormat::fmt(&direction, ref_direction_idx, f)
}

impl StepFormat for Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let transform = self.transform();
        let position_idx = idx + 1;
        let r0 = transform.column(0).into_owned().norm();
        let r1 = transform.column(1).into_owned().norm();
        if r0.near(&r1) {
            let r = FloatDisplay(r0);
            f.write_fmt(format_args!("#{idx} = CIRCLE('', #{position_idx}, {r});\n"))?;
        } else {
            let (r0, r1) = (FloatDisplay(r0), FloatDisplay(r1));
            f.write_fmt(format_args!(
                "#{idx} = ELLIPSE('', #{position_idx}, {r0}, {r1});\n"
            ))?;
        }
        write_conic_placement_3d(transform, position_idx, f)
    }
}
impl_const_step_length!(Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>, 5);

impl StepFormat for Processor<TrimmedCurve<UnitHyperbola<Point2>>, Matrix3> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let transform = *self.transform();
        let position_idx = idx + 1;
        let location_idx = idx + 2;
        let ref_direction_idx = idx + 3;
        let r0 = transform.column(0).into_owned().norm();
        let r1 = transform.column(1).into_owned().norm();
        let ref_direction_raw = VectorAsDirection(transform.column(0).into_owned().spatial() / r0);
        let ref_direction = StepDisplay::new(ref_direction_raw, ref_direction_idx);
        let location = StepDisplay::new(transform.column(2).into_owned().to_point(), location_idx);
        let (r0, r1) = (FloatDisplay(r0), FloatDisplay(r1));
        f.write_fmt(format_args!(
            "#{idx} = HYPERBOLA('', #{position_idx}, {r0}, {r1});
#{position_idx} = AXIS2_PLACEMENT_2D('', #{location_idx}, #{ref_direction_idx});
{location}{ref_direction}"
        ))
    }
}
impl_const_step_length!(Processor<TrimmedCurve<UnitHyperbola<Point2>>, Matrix3>, 4);

impl StepFormat for Processor<TrimmedCurve<UnitHyperbola<Point3>>, Matrix4> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let transform = self.transform();
        let position_idx = idx + 1;
        let r0 = transform.column(0).into_owned().norm();
        let r1 = transform.column(1).into_owned().norm();
        let (r0, r1) = (FloatDisplay(r0), FloatDisplay(r1));
        f.write_fmt(format_args!(
            "#{idx} = HYPERBOLA('', #{position_idx}, {r0}, {r1});\n"
        ))?;
        write_conic_placement_3d(transform, position_idx, f)
    }
}

impl_const_step_length!(Processor<TrimmedCurve<UnitHyperbola<Point3>>, Matrix4>, 5);

impl StepFormat for Processor<TrimmedCurve<UnitParabola<Point2>>, Matrix3> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let transform = *self.transform();
        let position_idx = idx + 1;
        let location_idx = idx + 2;
        let ref_direction_idx = idx + 3;
        let r0 = transform.column(0).into_owned().norm();
        let r1 = transform.column(1).into_owned().norm();
        let focal_dist = FloatDisplay(r1 * r1 / r0);
        let ref_direction_raw = VectorAsDirection(transform.column(0).into_owned().spatial() / r0);
        let ref_direction = StepDisplay::new(ref_direction_raw, ref_direction_idx);
        let location = StepDisplay::new(transform.column(2).into_owned().to_point(), location_idx);
        f.write_fmt(format_args!(
            "#{idx} = PARABOLA('', #{position_idx}, {focal_dist});
#{position_idx} = AXIS2_PLACEMENT_2D('', #{location_idx}, #{ref_direction_idx});
{location}{ref_direction}"
        ))
    }
}
impl_const_step_length!(Processor<TrimmedCurve<UnitParabola<Point2>>, Matrix3>, 4);

impl StepFormat for Processor<TrimmedCurve<UnitParabola<Point3>>, Matrix4> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let transform = self.transform();
        let position_idx = idx + 1;
        let r0 = transform.column(0).into_owned().norm();
        let r1 = transform.column(1).into_owned().norm();
        let focal_dist = FloatDisplay(r1 * r1 / r0);
        f.write_fmt(format_args!(
            "#{idx} = PARABOLA('', #{position_idx}, {focal_dist});\n"
        ))?;
        write_conic_placement_3d(transform, position_idx, f)
    }
}

impl_const_step_length!(Processor<TrimmedCurve<UnitParabola<Point3>>, Matrix4>, 5);

impl<C, M: One> StepCurve for Processor<C, M> {
    fn same_sense(&self) -> bool { self.orientation() }
}

impl<C, S0, S1> StepFormat for IntersectionCurve<C, S0, S1>
where
    C: StepLength + StepFormat,
    S0: StepLength + StepFormat,
    S1: StepLength + StepFormat,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let curve_idx = idx + 1;
        let surface0_idx = curve_idx + self.leader().step_length();
        let surface1_idx = surface0_idx + self.surface0().step_length();
        f.write_fmt(format_args!(
            "#{idx} = INTERSECTION_CURVE('', #{curve_idx}, (#{surface0_idx}, #{surface1_idx}), .CURVE_3D.);\n"
        ))?;
        self.leader().fmt(curve_idx, f)?;
        self.surface0().fmt(surface0_idx, f)?;
        self.surface1().fmt(surface1_idx, f)
    }
}

impl<C: StepLength, S0: StepLength, S1: StepLength> StepLength for IntersectionCurve<C, S0, S1> {
    fn step_length(&self) -> usize {
        1 + self.leader().step_length()
            + self.surface0().step_length()
            + self.surface1().step_length()
    }
}

impl<C, S0, S1, T0, T1> StepFormat for SurfaceCurve<C, S0, S1, T0, T1>
where
    C: StepLength + StepFormat,
    S0: StepLength + StepFormat,
    S1: StepLength + StepFormat,
    T0: StepLength + StepFormat,
    T1: StepLength + StepFormat,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let curve_idx = idx + 1;
        let boundary0_idx = curve_idx + self.leader().step_length();
        let boundary1_idx = boundary0_idx + self.boundary0().map_or(0, StepLength::step_length);
        let surface0_idx = boundary1_idx + self.boundary1().map_or(0, StepLength::step_length);
        let surface1_idx = surface0_idx + self.surface0().step_length();
        let assoc0 = self.boundary0().map_or_else(
            || format!("#{surface0_idx}"),
            |_| format!("#{boundary0_idx}"),
        );
        let assoc1 = self.boundary1().map_or_else(
            || format!("#{surface1_idx}"),
            |_| format!("#{boundary1_idx}"),
        );
        f.write_fmt(format_args!(
            "#{idx} = INTERSECTION_CURVE('', #{curve_idx}, ({assoc0}, {assoc1}), .CURVE_3D.);\n"
        ))?;
        self.leader().fmt(curve_idx, f)?;
        if let Some(boundary) = self.boundary0() {
            boundary.fmt(boundary0_idx, f)?;
        }
        if let Some(boundary) = self.boundary1() {
            boundary.fmt(boundary1_idx, f)?;
        }
        self.surface0().fmt(surface0_idx, f)?;
        self.surface1().fmt(surface1_idx, f)
    }
}

impl<C, S0, S1, T0, T1> StepLength for SurfaceCurve<C, S0, S1, T0, T1>
where
    C: StepLength,
    S0: StepLength,
    S1: StepLength,
    T0: StepLength,
    T1: StepLength,
{
    fn step_length(&self) -> usize {
        1 + self.leader().step_length()
            + self.boundary0().map_or(0, StepLength::step_length)
            + self.boundary1().map_or(0, StepLength::step_length)
            + self.surface0().step_length()
            + self.surface1().step_length()
    }
}

impl<C, S0, S1> ConstStepLength for IntersectionCurve<C, S0, S1>
where
    C: ConstStepLength,
    S0: ConstStepLength,
    S1: ConstStepLength,
{
    const LENGTH: usize = 1 + C::LENGTH + S0::LENGTH + S1::LENGTH;
}

impl<C: StepCurve, S0, S1> StepCurve for IntersectionCurve<C, S0, S1> {
    fn same_sense(&self) -> bool { self.leader().same_sense() }
}

// The parameter-space context specifies only a dimension and the fixed
// context identifier. The support and the actual chart coordinates belong
// to the PCURVE itself, not to this context.
pub(super) fn fmt_parameter_context(idx: usize, f: &mut Formatter<'_>) -> Result {
    f.write_fmt(format_args!(
        "#{idx} = (
    GEOMETRIC_REPRESENTATION_CONTEXT(2)
    PARAMETRIC_REPRESENTATION_CONTEXT()
    REPRESENTATION_CONTEXT('2D SPACE', '')
);
"
    ))
}

fn fmt_pcurve<C, S>(
    pcurve: &ParameterCurve<C, S>,
    idx: usize,
    support_reference: Option<usize>,
    context_reference: Option<usize>,
    f: &mut Formatter<'_>,
) -> Result
where
    C: StepFormat + Transformed<Matrix3> + Clone,
    S: StepFormat + StepLength,
{
    let surface_idx = support_reference.unwrap_or(idx + 1);
    let repr_idx =
        support_reference.map_or_else(|| surface_idx + pcurve.surface().step_length(), |_| idx + 1);
    let context_idx = context_reference.unwrap_or(repr_idx + 1);
    let curve_idx = repr_idx + 1 + usize::from(context_reference.is_none());
    let mapped = pcurve
        .surface()
        .parameter_space_transform()
        .map(|transform| pcurve.curve().transformed(transform));
    let curve = StepDisplay::new(mapped.as_ref().unwrap_or_else(|| pcurve.curve()), curve_idx);
    f.write_fmt(format_args!(
        "#{idx} = PCURVE('', #{surface_idx}, #{repr_idx});\n"
    ))?;
    if support_reference.is_none() {
        StepFormat::fmt(pcurve.surface(), surface_idx, f)?;
    }
    f.write_fmt(format_args!(
        "#{repr_idx} = DEFINITIONAL_REPRESENTATION('', (#{curve_idx}), #{context_idx});\n"
    ))?;
    if context_reference.is_none() {
        fmt_parameter_context(context_idx, f)?;
    }
    std::fmt::Display::fmt(&curve, f)
}

impl<C, S> StepFormat for ParameterCurve<C, S>
where
    C: StepFormat + Transformed<Matrix3> + Clone,
    S: StepFormat + StepLength,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        fmt_pcurve(self, idx, None, None, f)
    }
}

impl<C, S> StepPcurve for ParameterCurve<C, S>
where
    C: StepFormat + StepLength + Transformed<Matrix3> + Clone,
    S: StepFormat + StepLength,
{
    fn matches_surface_encoding(&self, encoding: &str) -> bool {
        canonical_step_encoding(self.surface()).is_some_and(|support| support == encoding)
    }

    fn step_length_with_shared_context(&self, shared_surface: bool) -> usize {
        2 + self.curve().step_length()
            + if shared_surface {
                0
            } else {
                self.surface().step_length()
            }
    }

    fn fmt_with_shared_context(
        &self,
        idx: usize,
        surface_idx: Option<usize>,
        context_idx: usize,
        f: &mut Formatter<'_>,
    ) -> Result {
        fmt_pcurve(self, idx, surface_idx, Some(context_idx), f)
    }
}

impl<C: StepLength, S: StepLength> StepLength for ParameterCurve<C, S> {
    fn step_length(&self) -> usize { 3 + self.curve().step_length() + self.surface().step_length() }
}

impl<C, S> ConstStepLength for ParameterCurve<C, S>
where
    C: ConstStepLength,
    S: ConstStepLength,
{
    const LENGTH: usize = 3 + C::LENGTH + S::LENGTH;
}

impl<C: StepCurve, S> StepCurve for ParameterCurve<C, S> {
    fn same_sense(&self) -> bool { self.curve().same_sense() }
}

impl StepFormat for ModelingConic2D {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Ellipse(x) => StepFormat::fmt(x, idx, f),
            Self::Hyperbola(x) => StepFormat::fmt(x, idx, f),
            Self::Parabola(x) => StepFormat::fmt(x, idx, f),
        }
    }
}

impl StepLength for ModelingConic2D {
    fn step_length(&self) -> usize {
        match self {
            Self::Ellipse(x) => x.step_length(),
            Self::Hyperbola(x) => x.step_length(),
            Self::Parabola(x) => x.step_length(),
        }
    }
}

impl StepCurve for ModelingConic2D {}

impl StepFormat for ModelingCurve2D {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Line(x) => StepFormat::fmt(x, idx, f),
            Self::Polyline(x) => StepFormat::fmt(x, idx, f),
            Self::Conic(x) => StepFormat::fmt(x, idx, f),
            Self::BsplineCurve(x) => StepFormat::fmt(x, idx, f),
            Self::NurbsCurve(x) => StepFormat::fmt(x, idx, f),
        }
    }
}

impl StepLength for ModelingCurve2D {
    fn step_length(&self) -> usize {
        match self {
            Self::Line(_) => Line::<Point2>::LENGTH,
            Self::Polyline(x) => x.step_length(),
            Self::Conic(x) => x.step_length(),
            Self::BsplineCurve(x) => x.step_length(),
            Self::NurbsCurve(x) => x.step_length(),
        }
    }
}

impl StepCurve for ModelingCurve2D {}

impl StepFormat for ModelingCurve {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Line(x) => StepFormat::fmt(x, idx, f),
            Self::BsplineCurve(x) => StepFormat::fmt(x, idx, f),
            Self::NurbsCurve(x) => StepFormat::fmt(x, idx, f),
            Self::ParameterCurve(x) => StepFormat::fmt(x, idx, f),
            Self::Intersection(x) => StepFormat::fmt(x, idx, f),
        }
    }
}

impl StepLength for ModelingCurve {
    fn step_length(&self) -> usize {
        match self {
            Self::Line(_) => Line::<Point3>::LENGTH,
            Self::BsplineCurve(x) => x.step_length(),
            Self::NurbsCurve(x) => x.step_length(),
            Self::ParameterCurve(x) => x.step_length(),
            Self::Intersection(x) => x.step_length(),
        }
    }
}

impl StepCurve for ModelingCurve {}

impl StepFormat for Plane {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let axis2_placement_idx = idx + 1;
        let location_idx = idx + 2;
        let z_axis_idx = idx + 3;
        let x_axis_idx = idx + 4;
        f.write_fmt(format_args!(
            "#{idx} = PLANE('', #{axis2_placement_idx});
#{axis2_placement_idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{z_axis_idx}, #{x_axis_idx});
{location}{z_axis}{x_axis}",
            location = StepDisplay::new(self.origin(), location_idx),
            z_axis = StepDisplay::new(VectorAsDirection(self.normal()), z_axis_idx),
            x_axis = StepDisplay::new(VectorAsDirection(self.axis_u().normalize()), x_axis_idx)
        ))
    }
}
impl_const_step_length!(Plane, 5);

impl StepSurface for Plane {}

/// A right-handed STEP elementary-surface placement from a certified
/// similarity. The center must undergo the full affine transformation,
/// not merely be added to the translation.
fn elementary_pose(center: Point3, transform: &Matrix4) -> Option<(Point3, Vector3, Vector3, f64)> {
    let scale = monstertruck_geometry::analytic_carrier::similarity_scale(*transform)?;
    let location = transform.apply_point(center);
    [location.x, location.y, location.z]
        .into_iter()
        .all(f64::is_finite)
        .then_some((
            location,
            transform.column(2).into_owned().spatial() / scale,
            transform.column(0).into_owned().spatial() / scale,
            scale,
        ))
}

fn sphere_pose(sphere: &Processor<Sphere, Matrix4>) -> Option<(Point3, Vector3, Vector3, f64)> {
    let (location, axis, reference, scale) =
        elementary_pose(sphere.entity().center(), sphere.transform())?;
    let radius = scale * sphere.entity().radius();
    (radius.is_finite() && radius > 0.0).then_some((location, axis, reference, radius))
}

impl StepFormat for Processor<Sphere, Matrix4> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let (location, axis, reference, radius) = sphere_pose(self).ok_or(std::fmt::Error)?;
        let position_idx = idx + 1;
        let location_idx = idx + 2;
        let axis_idx = idx + 3;
        let ref_direction_idx = idx + 4;
        let r = FloatDisplay(radius);
        f.write_fmt(format_args!(
            "#{idx} = SPHERICAL_SURFACE('', #{position_idx}, {r});
#{position_idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});
"
        ))?;
        StepFormat::fmt(&location, location_idx, f)?;
        StepFormat::fmt(&VectorAsDirection(axis), axis_idx, f)?;
        StepFormat::fmt(&VectorAsDirection(reference), ref_direction_idx, f)
    }
    fn parameter_space_transform(&self) -> Option<Matrix3> {
        sphere_pose(self)?;
        let sign = if self.transform().determinant() > 0.0 {
            1.0
        } else {
            -1.0
        };
        let longitude = Vector3::new(sign, 0.0, 0.0);
        let polar = Vector3::new(0.0, -1.0, 0.0);
        let translation = Vector3::new(0.0, std::f64::consts::FRAC_PI_2, 1.0);
        let (first, second) = if self.orientation() {
            (polar, longitude)
        } else {
            (longitude, polar)
        };
        Some(Matrix3::from_columns(&[first, second, translation]))
    }
}
impl_const_step_length!(Processor<Sphere, Matrix4>, 5);

impl StepFormat for Sphere {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        StepFormat::fmt(&Processor::new(*self), idx, f)
    }
}
impl_const_step_length!(Sphere, 5);
impl StepSurface for Sphere {}

fn torus_pose(torus: &Processor<Torus, Matrix4>) -> Option<(Point3, Vector3, Vector3, f64, f64)> {
    let (location, axis, reference, scale) =
        elementary_pose(torus.entity().center(), torus.transform())?;
    let major = scale * torus.entity().large_radius();
    let minor = scale * torus.entity().small_radius();
    ([major, minor]
        .iter()
        .all(|radius| radius.is_finite() && *radius > 0.0))
    .then_some((location, axis, reference, major, minor))
}

impl StepFormat for Processor<Torus, Matrix4> {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let (location, axis, reference, major, minor) = torus_pose(self).ok_or(std::fmt::Error)?;
        let position_idx = idx + 1;
        let location_idx = idx + 2;
        let axis_idx = idx + 3;
        let ref_direction_idx = idx + 4;
        let greater = FloatDisplay(major);
        let lesser = FloatDisplay(minor);
        f.write_fmt(format_args!(
            "#{idx} = TOROIDAL_SURFACE('', #{position_idx}, {greater}, {lesser});
#{position_idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});
",
        ))?;
        StepFormat::fmt(&location, location_idx, f)?;
        StepFormat::fmt(&VectorAsDirection(axis), axis_idx, f)?;
        StepFormat::fmt(&VectorAsDirection(reference), ref_direction_idx, f)
    }
    fn parameter_space_transform(&self) -> Option<Matrix3> {
        torus_pose(self)?;
        let sign = if self.transform().determinant() > 0.0 {
            1.0
        } else {
            -1.0
        };
        let longitude = Vector3::new(sign, 0.0, 0.0);
        let meridian = Vector3::new(0.0, 1.0, 0.0);
        let (first, second) = if self.orientation() {
            (longitude, meridian)
        } else {
            (meridian, longitude)
        };
        Some(Matrix3::from_columns(&[first, second, Vector3::z()]))
    }
}
impl_const_step_length!(Processor<Torus, Matrix4>, 5);

impl StepSurface for Processor<Torus, Matrix4> {
    fn same_sense(&self) -> bool { self.orientation() == (self.transform().determinant() > 0.0) }
}

impl StepFormat for Torus {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        StepFormat::fmt(&Processor::new(*self), idx, f)
    }
}
impl_const_step_length!(Torus, 5);
impl StepSurface for Torus {}

impl<P> StepFormat for BsplineSurface<P>
where P: Copy + StepFormat
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let control_points = self.control_points();
        let control_points_instances = self
            .control_points()
            .iter()
            .flatten()
            .enumerate()
            .map(|(i, p)| StepDisplay::new(*p, idx + i + 1))
            .collect::<Vec<_>>();
        let mut counter = 0;
        let control_points_list = control_points
            .iter()
            .map(|slice| {
                counter += slice.len();
                IndexSliceDisplay(idx + counter - slice.len() + 1..=idx + counter)
            })
            .collect::<Vec<_>>();
        let (horizontal_knots, horizontal_multiplicities) = self.knot_vector_u().to_single_multi();
        let (vertical_knots, vertical_multiplicities) = self.knot_vector_v().to_single_multi();
        f.write_fmt(format_args!(
            "#{idx} = B_SPLINE_SURFACE_WITH_KNOTS('', {u_degree}, {v_degree}, {control_points_list}, .UNSPECIFIED., .U., .U., .U., \
{u_multiplicities}, {v_multiplicities}, {u_knots}, {v_knots}, .UNSPECIFIED.);\n{control_points_instances}",
            u_degree = self.udegree(),
            v_degree = self.vdegree(),
            control_points_list = SliceDisplay(&control_points_list),
            u_multiplicities = SliceDisplay(&horizontal_multiplicities),
            v_multiplicities = SliceDisplay(&vertical_multiplicities),
            u_knots = SliceDisplay(&horizontal_knots),
            v_knots = SliceDisplay(&vertical_knots),
            control_points_instances = SliceDisplay(&control_points_instances),
        ))
    }
}

impl<P> StepLength for BsplineSurface<P> {
    fn step_length(&self) -> usize { 1 + self.control_points().iter().map(Vec::len).sum::<usize>() }
}
impl<P> StepSurface for BsplineSurface<P> {}

impl<V> StepFormat for NurbsSurface<V>
where
    V: Homogeneous<Scalar = f64>,
    V::Point: Copy + StepFormat,
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let control_points_instances = self
            .control_points()
            .iter()
            .flatten()
            .enumerate()
            .map(|(i, v)| StepDisplay::new(v.to_point(), idx + i + 1))
            .collect::<Vec<_>>();
        let mut counter = 0;
        let control_points_list = self
            .control_points()
            .iter()
            .map(|slice| {
                counter += slice.len();
                IndexSliceDisplay(idx + counter - slice.len() + 1..=idx + counter)
            })
            .collect::<Vec<_>>();
        let weights = self
            .control_points()
            .iter()
            .map(|slice| slice.iter().map(|v| v.weight()).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let weights = weights
            .iter()
            .map(|slice| SliceDisplay(slice))
            .collect::<Vec<_>>();
        let (horizontal_knots, horizontal_multiplicities) = self.knot_vector_u().to_single_multi();
        let (vertical_knots, vertical_multiplicities) = self.knot_vector_v().to_single_multi();
        f.write_fmt(format_args!(
            "#{idx} = (
    BOUNDED_SURFACE()
    B_SPLINE_SURFACE({u_degree}, {v_degree}, {control_points_list}, .UNSPECIFIED., .U., .U., .U.)
    B_SPLINE_SURFACE_WITH_KNOTS({u_multiplicities}, {v_multiplicities}, {u_knots}, {v_knots}, .UNSPECIFIED.)
    GEOMETRIC_REPRESENTATION_ITEM()
    RATIONAL_B_SPLINE_SURFACE({weights})
    REPRESENTATION_ITEM('')
    SURFACE()
);\n{control_points_instances}",
            u_degree = self.udegree(),
            v_degree = self.vdegree(),
            control_points_list = SliceDisplay(&control_points_list),
            u_multiplicities = SliceDisplay(&horizontal_multiplicities),
            v_multiplicities = SliceDisplay(&vertical_multiplicities),
            u_knots = SliceDisplay(&horizontal_knots),
            v_knots = SliceDisplay(&vertical_knots),
            control_points_instances = SliceDisplay(&control_points_instances),
            weights = SliceDisplay(&weights),
        ))
    }
}

impl<V> StepLength for NurbsSurface<V> {
    fn step_length(&self) -> usize { 1 + self.control_points().iter().map(Vec::len).sum::<usize>() }
}
impl<V> StepSurface for NurbsSurface<V> {}

impl<C> StepFormat for ExtrusionSurface<C, Vector3>
where C: StepLength + StepFormat
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let curve = self.entity_curve();
        let curve_idx = idx + 1;
        let vector_idx = idx + 1 + curve.step_length();
        let vector = self.extruding_vector();
        f.write_fmt(format_args!(
            "#{idx} = SURFACE_OF_LINEAR_EXTRUSION('', #{curve_idx}, #{vector_idx});\n{}{}",
            StepDisplay::new(curve, curve_idx),
            StepDisplay::new(vector, vector_idx),
        ))
    }
}
impl<C: StepLength> StepLength for ExtrusionSurface<C, Vector3> {
    fn step_length(&self) -> usize { 1 + self.entity_curve().step_length() + Vector3::LENGTH }
}
impl<C: ConstStepLength> ConstStepLength for ExtrusionSurface<C, Vector3> {
    const LENGTH: usize = 1 + C::LENGTH + Vector3::LENGTH;
}
impl<C> StepSurface for ExtrusionSurface<C, Vector3> {}

impl<C, T: One> StepSurface for Processor<ExtrusionSurface<C, Vector3>, T> {
    fn same_sense(&self) -> bool { self.orientation() }
}

impl<C> StepFormat for RevolutionSurface<C>
where C: StepLength + StepFormat
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let curve = self.entity_curve();
        let curve_idx = idx + 1;
        let axis_idx = curve_idx + curve.step_length();
        let location_idx = axis_idx + 1;
        let dir_idx = location_idx + 1;
        f.write_fmt(format_args!(
            "#{idx} = SURFACE_OF_REVOLUTION('', #{curve_idx}, #{axis_idx});
{curve}#{axis_idx} = AXIS1_PLACEMENT('', #{location_idx}, #{dir_idx});\n{location}{dir}",
            curve = StepDisplay::new(curve, curve_idx),
            location = StepDisplay::new(self.origin(), location_idx),
            dir = StepDisplay::new(VectorAsDirection(self.axis()), dir_idx),
        ))
    }
}

impl<C: StepLength> StepLength for RevolutionSurface<C> {
    fn step_length(&self) -> usize { 4 + self.entity_curve().step_length() }
}

impl<C: ConstStepLength> ConstStepLength for RevolutionSurface<C> {
    const LENGTH: usize = 4 + C::LENGTH;
}

impl<C> StepSurface for RevolutionSurface<C> {
    fn same_sense(&self) -> bool { false }
}

impl<C> StepFormat for Processor<RevolutionSurface<C>, Matrix4>
where C: StepLength + Transformed<Matrix4> + StepFormat
{
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let surface = self.entity();
        let transform = self.transform();
        let Some((k, a)) = monstertruck_core::numerical::affine_qr(*transform) else {
            f.write_str("Transform is not regular")?;
            return ERR;
        };
        if !a[(0, 0)].near(&a[(1, 1)]) || !a[(1, 1)].near(&a[(2, 2)]) {
            f.write_str("Transform contains non-uniform scale.")?;
            return ERR;
        }
        let curve = surface.entity_curve().transformed(*transform);
        let axis = k.apply_vector(surface.axis());
        let origin = transform.apply_point(surface.origin());
        let surface = RevolutionSurface::by_revolution(curve, origin, axis);
        StepFormat::fmt(&surface, idx, f)
    }
}
impl<C: StepLength> StepLength for Processor<RevolutionSurface<C>, Matrix4> {
    fn step_length(&self) -> usize { self.entity().step_length() }
}

impl<C, T: One> StepSurface for Processor<RevolutionSurface<C>, T> {
    fn same_sense(&self) -> bool { !self.orientation() }
}

/// Placement of a STEP `CYLINDRICAL_SURFACE`, expressed in world coordinates.
struct CylindricalSurfaceData {
    location: Point3,
    axis: Vector3,
    ref_direction: Vector3,
    radius: f64,
    /// Native parameter chart → right-handed STEP (azimuth, axial length).
    uv: Matrix3,
}

impl StepFormat for CylindricalSurfaceData {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        let position_idx = idx + 1;
        let location_idx = idx + 2;
        let axis_idx = idx + 3;
        let ref_direction_idx = idx + 4;
        let radius = FloatDisplay(self.radius);
        f.write_fmt(format_args!(
            "#{idx} = CYLINDRICAL_SURFACE('', #{position_idx}, {radius});
#{position_idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});\n"
        ))?;
        StepFormat::fmt(&self.location, location_idx, f)?;
        StepFormat::fmt(&VectorAsDirection(self.axis), axis_idx, f)?;
        StepFormat::fmt(&VectorAsDirection(self.ref_direction), ref_direction_idx, f)
    }
}
impl_const_step_length!(CylindricalSurfaceData, 5);

/// Recognizes a right circular cylinder expressed as a surface of revolution of
/// a straight profile line parallel to the revolution axis, returning its STEP
/// `CYLINDRICAL_SURFACE` placement in world coordinates. Returns `None` for any
/// other revolution (cones, spheres, general profiles) or a transform that does
/// not preserve the circular cross-section.
fn cylindrical_surface_from_revolution(
    surface: &Processor<RevolutionSurface<ModelingCurve>, Matrix4>,
) -> Option<CylindricalSurfaceData> {
    let ModelingCurve::Line(profile) = surface.entity().entity_curve() else {
        return None;
    };
    let transform = surface.transform();
    let scale = monstertruck_geometry::analytic_carrier::similarity_scale(*transform)?;
    let revolution = surface.entity();
    let axis = revolution.axis().normalize();
    let direction = profile.1 - profile.0;
    // The profile must be parallel to the axis so the radius stays constant.
    if direction.so_small() || !direction.normalize().cross(&axis).so_small() {
        return None;
    }
    let offset = profile.0 - revolution.origin();
    let radial = offset - axis * offset.dot(&axis);
    if radial.so_small() {
        return None;
    }
    let foot = revolution.origin() + axis * offset.dot(&axis);
    let radial = transform.apply_vector(radial);
    let angular_sign = transform.determinant().signum();
    let axial = Vector3::new(0.0, direction.dot(&axis) * scale, 0.0);
    let angular = Vector3::new(angular_sign, 0.0, 0.0);
    let (first, second) = if surface.orientation() {
        (axial, angular)
    } else {
        (angular, axial)
    };
    let uv = Matrix3::from_columns(&[first, second, Vector3::z()]);
    Some(CylindricalSurfaceData {
        location: transform.apply_point(foot),
        axis: transform.apply_vector(axis).normalize(),
        ref_direction: radial.normalize(),
        radius: radial.norm(),
        uv,
    })
}

impl StepFormat for ModelingSurface {
    fn fmt(&self, idx: usize, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Plane(x) => StepFormat::fmt(x, idx, f),
            Self::BsplineSurface(x) => StepFormat::fmt(x, idx, f),
            Self::NurbsSurface(x) => StepFormat::fmt(x, idx, f),
            Self::Revolution(x) => {
                if let Some(cylinder) = cylindrical_surface_from_revolution(x) {
                    StepFormat::fmt(&cylinder, idx, f)
                } else if let Some(bsp) = x.try_into_homogeneous_bspline_surface() {
                    let nurbs = NurbsSurface::new(bsp);
                    StepFormat::fmt(&nurbs, idx, f)
                } else {
                    StepFormat::fmt(x, idx, f)
                }
            }
            Self::Tspline(t_mesh) => {
                let bsp = t_mesh.to_bspline_surface(4).map_err(|_| std::fmt::Error)?;
                StepFormat::fmt(&bsp, idx, f)
            }
            Self::Sphere(sphere) if sphere_pose(sphere).is_some() => {
                StepFormat::fmt(sphere, idx, f)
            }
            Self::Torus(torus) if torus_pose(torus).is_some() => StepFormat::fmt(torus, idx, f),
            // Non-similarity transforms require the full rational net.
            Self::Sphere(_) | Self::Torus(_) => self
                .try_into_homogeneous_bspline_surface()
                .map_or(Err(std::fmt::Error), |bsp| {
                    StepFormat::fmt(&NurbsSurface::new(bsp), idx, f)
                }),
        }
    }
    fn parameter_space_transform(&self) -> Option<Matrix3> {
        match self {
            Self::Revolution(revolution) => {
                cylindrical_surface_from_revolution(revolution).map(|cylinder| cylinder.uv)
            }
            Self::Sphere(sphere) => sphere.parameter_space_transform(),
            Self::Torus(torus) => torus.parameter_space_transform(),
            _ => None,
        }
    }
}

impl StepLength for ModelingSurface {
    fn step_length(&self) -> usize {
        match self {
            Self::Plane(_) => Plane::LENGTH,
            Self::BsplineSurface(x) => x.step_length(),
            Self::NurbsSurface(x) => x.step_length(),
            Self::Revolution(x) => {
                if cylindrical_surface_from_revolution(x).is_some() {
                    CylindricalSurfaceData::LENGTH
                } else if let Some(bsp) = x.try_into_homogeneous_bspline_surface() {
                    NurbsSurface::new(bsp).step_length()
                } else {
                    x.entity().step_length()
                }
            }
            Self::Tspline(t_mesh) => t_mesh
                .to_bspline_surface(4)
                .map_or(0, |bsp| bsp.step_length()),
            Self::Sphere(sphere) if sphere_pose(sphere).is_some() => sphere.step_length(),
            Self::Torus(torus) if torus_pose(torus).is_some() => torus.step_length(),
            Self::Sphere(_) | Self::Torus(_) => self
                .try_into_homogeneous_bspline_surface()
                .map_or(0, |bsp| NurbsSurface::new(bsp).step_length()),
        }
    }
}

impl StepSurface for ModelingSurface {
    fn same_sense(&self) -> bool {
        match self {
            // Cylinder placement always defines a right-handed outward
            // STEP chart. Flipping its axis alone cannot invert its normal;
            // the oriented chart determinant must propagate to ADVANCED_FACE.
            Self::Revolution(revolution) => cylindrical_surface_from_revolution(revolution)
                .is_none_or(|cylinder| cylinder.uv.determinant() > 0.0),
            Self::Sphere(sphere) if sphere_pose(sphere).is_some() => {
                sphere.orientation() == (sphere.transform().determinant() > 0.0)
            }
            Self::Torus(torus) if torus_pose(torus).is_some() => StepSurface::same_sense(torus),
            _ => true,
        }
    }
}
