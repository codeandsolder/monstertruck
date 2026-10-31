//! Surface entities: elementary surfaces, B-spline/NURBS surfaces and
//! swept surfaces.

use super::{
    Axis1Placement, Axis1PlacementHolder, Axis2Placement3d, Axis2Placement3dHolder, BsplineSurface,
    CartesianPoint, CartesianPointHolder, Curve3D, CurveAny, CurveAnyHolder, Deserialize,
    ElementarySurface, EuclideanSpace, Holder, InnerSpace, Invertible, KERNEL_EPSILON, KnotType,
    KnotVector, Logical, Matrix4, NurbsSurface, PI, Point3, Result, RevolutionSurface, Serialize,
    Sphere, StepConvertingError, StepExtrusionSurface, StepRevolutionSurface, Surface,
    SweepSurface, Table, Torus, Transformed, Vector, Vector3, Vector4, VectorHolder, geom,
    quasi_uniform_knots, step_geometry, uniform_knots,
};

/// `surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(generate_deserialize)]
pub enum SurfaceAny {
    #[holder(use_place_holder)]
    Elementary(Box<ElementarySurfaceAny>),
    #[holder(use_place_holder)]
    Bspline(Box<BsplineSurfaceAny>),
    #[holder(use_place_holder)]
    Swept(Box<SweptSurfaceAny>),
}

impl TryFrom<&SurfaceAny> for Surface {
    type Error = StepConvertingError;
    fn try_from(x: &SurfaceAny) -> Result<Self, Self::Error> {
        use SurfaceAny::{Bspline, Elementary, Swept};
        Ok(match x {
            Elementary(x) => Self::Elementary(x.as_ref().try_into()?),
            Bspline(x) => x.as_ref().try_into()?,
            Swept(x) => Self::Sweep(x.as_ref().try_into()?),
        })
    }
}

/// `elementary_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(generate_deserialize)]
pub enum ElementarySurfaceAny {
    #[holder(use_place_holder)]
    Plane(Box<Plane>),
    #[holder(use_place_holder)]
    SphericalSurface(Box<SphericalSurface>),
    #[holder(use_place_holder)]
    CylindricalSurface(Box<CylindricalSurface>),
    #[holder(use_place_holder)]
    ToroidalSurface(Box<ToroidalSurface>),
    /// STEP's sheet-selecting degenerate torus subtype. Present so the record
    /// deserializes into a named class instead of `Table::dummy`; true spindle
    /// instances remain typed refusals until the selected sheet is represented.
    #[holder(use_place_holder)]
    DegenerateToroidalSurface(Box<DegenerateToroidalSurface>),
    #[holder(use_place_holder)]
    ConicalSurface(Box<ConicalSurface>),
}

/// Fallible because STEP input can carry non-finite/non-positive elementary
/// surface parameters and because `DEGENERATE_TOROIDAL_SURFACE` has sheet
/// semantics that a plain full torus does not encode. Invalid scalar domains and
/// unsupported sheet selection therefore refuse typed before reaching geometry
/// constructors.
impl TryFrom<&ElementarySurfaceAny> for ElementarySurface {
    type Error = StepConvertingError;
    fn try_from(value: &ElementarySurfaceAny) -> Result<Self, Self::Error> {
        use ElementarySurfaceAny::{
            ConicalSurface, CylindricalSurface, DegenerateToroidalSurface, Plane, SphericalSurface,
            ToroidalSurface,
        };
        Ok(match value {
            Plane(x) => Self::Plane(x.as_ref().into()),
            SphericalSurface(x) => Self::Sphere(x.as_ref().into()),
            CylindricalSurface(x) => Self::CylindricalSurface(x.as_ref().into()),
            ToroidalSurface(x) => Self::ToroidalSurface(x.as_ref().try_into()?),
            DegenerateToroidalSurface(x) => Self::ToroidalSurface(x.as_ref().try_into()?),
            ConicalSurface(x) => Self::ConicalSurface(x.as_ref().into()),
        })
    }
}

/// `plane`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = plane)]
#[holder(generate_deserialize)]
pub struct Plane {
    label: String,
    #[holder(use_place_holder)]
    position: Axis2Placement3d,
}

impl From<&Plane> for geom::Plane {
    fn from(plane: &Plane) -> Self {
        let mat = Matrix4::from(&plane.position);
        let o = Point3::from_homogeneous(mat[3]);
        let p = o + mat[0].truncate();
        let q = o + mat[1].truncate();
        Self::new(o, p, q)
    }
}

/// `spherical_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = spherical_surface)]
#[holder(generate_deserialize)]
pub struct SphericalSurface {
    label: String,
    #[holder(use_place_holder)]
    position: Axis2Placement3d,
    radius: f64,
}

impl From<&SphericalSurface> for step_geometry::SphericalSurface {
    fn from(ss: &SphericalSurface) -> Self {
        let mat = Matrix4::from(&ss.position);
        let sphere = Sphere(geom::Sphere::new(Point3::origin(), ss.radius));
        Self::new(sphere).transformed(mat)
    }
}

/// `cylindrical_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = cylindrical_surface)]
#[holder(generate_deserialize)]
pub struct CylindricalSurface {
    label: String,
    #[holder(use_place_holder)]
    position: Axis2Placement3d,
    radius: f64,
}

impl From<&CylindricalSurface> for step_geometry::CylindricalSurface {
    fn from(cs: &CylindricalSurface) -> Self {
        let mat = Matrix4::from(&cs.position);
        let x = mat[0].truncate();
        let z = mat[2].truncate();
        let center = Point3::from_homogeneous(mat[3]);
        let radius = cs.radius;
        let p = center + x * radius;
        let mut res = Self::new(RevolutionSurface::by_revolution(
            geom::Line(p, p + z),
            center,
            z,
        ));
        res.invert();
        res
    }
}

/// `toroidal_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = toroidal_surface)]
#[holder(generate_deserialize)]
pub struct ToroidalSurface {
    label: String,
    #[holder(use_place_holder)]
    position: Axis2Placement3d,
    major_radius: f64,
    minor_radius: f64,
}

/// Validate the scalar domain shared by the STEP torus spellings.
///
/// A positive-major spindle (`0 < R < r`) is a perfectly well-defined
/// parametric surface. Monstertruck used to refuse it only because `Torus` could
/// not invert the negative signed-radial branch; the analytic inverse is now
/// branch-complete, so ordinary `TOROIDAL_SURFACE` records may use the full
/// positive `(R, r)` domain. A negative major radius remains out of schema and
/// would require an explicit UV reparameterisation of any source pcurves.
fn validate_torus_radii(
    entity: &str,
    major_radius: f64,
    minor_radius: f64,
) -> Result<(), StepConvertingError> {
    if !major_radius.is_finite() || !minor_radius.is_finite() {
        return Err(format!(
            "{entity} refused: non-finite radii (major_radius {major_radius}, minor_radius {minor_radius})."
        )
        .into());
    }
    if minor_radius <= 0.0 {
        return Err(format!(
            "{entity} refused: non-positive minor_radius {minor_radius} (major_radius {major_radius}); ISO 10303-42 declares it a positive_length_measure."
        )
        .into());
    }
    if major_radius <= 0.0 {
        return Err(format!(
            "{entity} refused: non-positive major_radius {major_radius} with minor_radius {minor_radius}; ISO 10303-42 declares major_radius a positive_length_measure. Supporting this exporter extension safely requires carrying the equivalent (u + pi, pi - v) UV reparameterisation into face-local trims."
        )
        .into());
    }
    Ok(())
}

impl TryFrom<&ToroidalSurface> for step_geometry::ToroidalSurface {
    type Error = StepConvertingError;
    fn try_from(
        ToroidalSurface {
            position,
            major_radius,
            minor_radius,
            ..
        }: &ToroidalSurface,
    ) -> Result<Self, Self::Error> {
        validate_torus_radii("TOROIDAL_SURFACE", *major_radius, *minor_radius)?;
        let mat = Matrix4::from(position);
        let torus = Torus::new(Point3::origin(), *major_radius, *minor_radius)?;
        Ok(Self::new(torus).transformed(mat))
    }
}

/// `degenerate_toroidal_surface`
///
/// ISO 10303-42 subtype of `toroidal_surface`:
/// `(name, position, major_radius, minor_radius, select_outer)`. The fifth
/// attribute is a BOOLEAN -- which SHEET of the self-intersecting torus the face
/// lies on -- and is nothing like `toroidal_surface`'s attribute list, which
/// simply stops at `minor_radius`. The subtype's WHERE rule is
/// `major_radius < minor_radius`.
///
/// The entity exists in the schema, so it deserializes here rather than falling
/// into `Table::dummy`: a record that lands in `dummy` can only ever be refused
/// as a lookup miss (`Lookup failed for #145`), which names neither the class nor
/// the reason. Conversion keeps the subtype distinct so its selected-sheet
/// semantics cannot silently become a full immersed torus.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = degenerate_toroidal_surface)]
#[holder(generate_deserialize)]
pub struct DegenerateToroidalSurface {
    label: String,
    #[holder(use_place_holder)]
    position: Axis2Placement3d,
    major_radius: f64,
    minor_radius: f64,
    select_outer: bool,
}

impl TryFrom<&DegenerateToroidalSurface> for step_geometry::ToroidalSurface {
    type Error = StepConvertingError;
    fn try_from(
        DegenerateToroidalSurface {
            position,
            major_radius,
            minor_radius,
            select_outer,
            ..
        }: &DegenerateToroidalSurface,
    ) -> Result<Self, Self::Error> {
        validate_torus_radii("DEGENERATE_TOROIDAL_SURFACE", *major_radius, *minor_radius)?;
        // Unlike plain TOROIDAL_SURFACE, the STEP subtype carries `select_outer`:
        // in the true spindle regime the entity denotes ONE sheet, not the whole
        // immersed torus. Until that sheet domain is represented explicitly,
        // accepting it as a full Torus would add geometry the file did not name.
        if *minor_radius - *major_radius > KERNEL_EPSILON * (*major_radius + *minor_radius) {
            return Err(format!(
                "DEGENERATE_TOROIDAL_SURFACE refused: spindle sheet selection (major_radius {major_radius} < minor_radius {minor_radius}, select_outer = {select_outer}) needs an explicit sheet domain; the full analytic spindle surface is representable but would not preserve this entity's semantics."
            )
            .into());
        }
        // At the horn limit the inner sheet collapses to the axis point. A ring
        // surface mis-spelled as degenerate is likewise only meaningful here for
        // the outer selection.
        if !select_outer {
            return Err(format!(
                "DEGENERATE_TOROIDAL_SURFACE refused: select_outer = .F. selects the inner \
                 (apex) sheet, which at major_radius {major_radius} / minor_radius \
                 {minor_radius} has collapsed onto the axis of revolution and is not a \
                 representable surface."
            )
            .into());
        }
        let mat = Matrix4::from(position);
        let torus = Torus::new(Point3::origin(), *major_radius, *minor_radius)?;
        Ok(Self::new(torus).transformed(mat))
    }
}

/// `conical_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = conical_surface)]
#[holder(generate_deserialize)]
pub struct ConicalSurface {
    label: String,
    #[holder(use_place_holder)]
    position: Axis2Placement3d,
    radius: f64,
    semi_angle: f64,
}

fn normalize_conical_semi_angle(semi_angle: f64) -> f64 {
    if semi_angle.abs() > PI {
        semi_angle.to_radians()
    } else {
        semi_angle
    }
}

impl From<&ConicalSurface> for step_geometry::ConicalSurface {
    fn from(
        ConicalSurface {
            position,
            radius,
            semi_angle,
            ..
        }: &ConicalSurface,
    ) -> Self {
        let mat = Matrix4::from(position);
        let p = Point3::new(*radius, 0.0, 0.0);
        let semi_angle = normalize_conical_semi_angle(*semi_angle);
        let v = Vector3::new(f64::tan(semi_angle), 0.0, 1.0);
        let rev = RevolutionSurface::by_revolution(
            geom::Line(p, p + v),
            Point3::origin(),
            Vector3::unit_z(),
        );
        let mut processor = Self::new(rev);
        processor.transform_by(mat);
        processor.invert();
        processor
    }
}

/// `b_spline_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(generate_deserialize)]
pub enum BsplineSurfaceAny {
    #[holder(use_place_holder)]
    NonRationalBsplineSurface(NonRationalBsplineSurface),
    #[holder(use_place_holder)]
    RationalBsplineSurface(RationalBsplineSurface),
}

impl TryFrom<&BsplineSurfaceAny> for Surface {
    type Error = StepConvertingError;
    fn try_from(value: &BsplineSurfaceAny) -> Result<Self, Self::Error> {
        use BsplineSurfaceAny::{NonRationalBsplineSurface, RationalBsplineSurface};
        Ok(match value {
            NonRationalBsplineSurface(bsp) => Self::Bspline(bsp.try_into()?),
            RationalBsplineSurface(bsp) => Self::Nurbs(bsp.try_into()?),
        })
    }
}

/// `b_spline_surface_form`
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum BsplineSurfaceForm {
    PlaneSurf,
    CylindricalSurf,
    ConicalSurf,
    SphericalSurf,
    ToroidalSurf,
    SurfOfRevolution,
    RuledSurf,
    GeneralisedCone,
    QuadricSurf,
    SurfOfLinearExtrusion,
    Unspecified,
}

/// `b_spline_surface_with_knots`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = b_spline_surface_with_knots)]
#[holder(generate_deserialize)]
pub struct BsplineSurfaceWithKnots {
    label: String,
    u_degree: i64,
    v_degree: i64,
    #[holder(use_place_holder)]
    control_points_list: Vec<Vec<CartesianPoint>>,
    surface_form: BsplineSurfaceForm,
    u_closed: Logical,
    v_closed: Logical,
    self_intersect: Logical,
    u_multiplicities: Vec<i64>,
    v_multiplicities: Vec<i64>,
    u_knots: Vec<f64>,
    v_knots: Vec<f64>,
    knot_spec: KnotType,
}

impl TryFrom<&BsplineSurfaceWithKnots> for BsplineSurface<Point3> {
    type Error = StepConvertingError;
    fn try_from(surface: &BsplineSurfaceWithKnots) -> Result<Self, StepConvertingError> {
        let horizontal_knots = surface.u_knots.clone();
        let horizontal_multiplicities = surface
            .u_multiplicities
            .iter()
            .map(|n| usize::try_from(*n))
            .collect::<Result<Vec<_>, _>>()?;
        let horizontal_knots =
            KnotVector::from_single_multi(horizontal_knots, horizontal_multiplicities)?;
        let vertical_knots = surface.v_knots.clone();
        let vertical_multiplicities = surface
            .v_multiplicities
            .iter()
            .map(|n| usize::try_from(*n))
            .collect::<Result<Vec<_>, _>>()?;
        let vertical_knots =
            KnotVector::from_single_multi(vertical_knots, vertical_multiplicities)?;
        let ctrls = surface
            .control_points_list
            .iter()
            .map(|vec| vec.iter().map(Point3::from).collect())
            .collect();
        Ok(Self::new((horizontal_knots, vertical_knots), ctrls)?)
    }
}

/// `uniform_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = uniform_surface)]
#[holder(generate_deserialize)]
pub struct UniformSurface {
    label: String,
    u_degree: i64,
    v_degree: i64,
    #[holder(use_place_holder)]
    control_points_list: Vec<Vec<CartesianPoint>>,
    surface_form: BsplineSurfaceForm,
    u_closed: Logical,
    v_closed: Logical,
    self_intersect: Logical,
}

impl TryFrom<&UniformSurface> for BsplineSurface<Point3> {
    type Error = StepConvertingError;
    fn try_from(surface: &UniformSurface) -> Result<Self, StepConvertingError> {
        let horizontal_knots = uniform_knots(
            surface.control_points_list.len(),
            usize::try_from(surface.u_degree)?,
        )?;
        let first = surface
            .control_points_list
            .first()
            .ok_or("control points list is empty.")?;
        let vertical_knots = uniform_knots(first.len(), usize::try_from(surface.v_degree)?)?;
        let ctrls = surface
            .control_points_list
            .iter()
            .map(|vec| vec.iter().map(Point3::from).collect())
            .collect();
        Ok(Self::new((horizontal_knots, vertical_knots), ctrls)?)
    }
}

/// `quasi_uniform_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = quasi_uniform_surface)]
#[holder(generate_deserialize)]
pub struct QuasiUniformSurface {
    label: String,
    u_degree: i64,
    v_degree: i64,
    #[holder(use_place_holder)]
    control_points_list: Vec<Vec<CartesianPoint>>,
    surface_form: BsplineSurfaceForm,
    u_closed: Logical,
    v_closed: Logical,
    self_intersect: Logical,
}

impl TryFrom<&QuasiUniformSurface> for BsplineSurface<Point3> {
    type Error = StepConvertingError;
    fn try_from(surface: &QuasiUniformSurface) -> Result<Self, StepConvertingError> {
        let horizontal_knots = quasi_uniform_knots(
            surface.control_points_list.len(),
            usize::try_from(surface.u_degree)?,
        )?;
        let first = surface
            .control_points_list
            .first()
            .ok_or("control points list is empty.")?;
        let vertical_knots = quasi_uniform_knots(first.len(), usize::try_from(surface.v_degree)?)?;
        let ctrls = surface
            .control_points_list
            .iter()
            .map(|vec| vec.iter().map(Point3::from).collect())
            .collect();
        Ok(Self::new((horizontal_knots, vertical_knots), ctrls)?)
    }
}

/// `bezier_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = bezier_surface)]
#[holder(generate_deserialize)]
pub struct BezierSurface {
    label: String,
    u_degree: i64,
    v_degree: i64,
    #[holder(use_place_holder)]
    control_points_list: Vec<Vec<CartesianPoint>>,
    surface_form: BsplineSurfaceForm,
    u_closed: Logical,
    v_closed: Logical,
    self_intersect: Logical,
}

impl TryFrom<&BezierSurface> for BsplineSurface<Point3> {
    type Error = StepConvertingError;

    fn try_from(value: &BezierSurface) -> Result<Self, Self::Error> {
        let horizontal_knots = KnotVector::bezier_knot(usize::try_from(value.u_degree)?);
        let vertical_knots = KnotVector::bezier_knot(usize::try_from(value.v_degree)?);
        let ctrls = value
            .control_points_list
            .iter()
            .map(|vec| vec.iter().map(Point3::from).collect())
            .collect();
        Self::new((horizontal_knots, vertical_knots), ctrls).map_err(Into::into)
    }
}

/// Entity that does not exist in AP042.
/// Surface before rationalization of [`RationalBsplineSurface`] defined by a complex entity
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(generate_deserialize)]
pub enum NonRationalBsplineSurface {
    #[holder(use_place_holder)]
    BsplineSurfaceWithKnots(Box<BsplineSurfaceWithKnots>),
    #[holder(use_place_holder)]
    UniformSurface(Box<UniformSurface>),
    #[holder(use_place_holder)]
    QuasiUniformSurface(Box<QuasiUniformSurface>),
    #[holder(use_place_holder)]
    BezierSurface(Box<BezierSurface>),
}

impl TryFrom<&NonRationalBsplineSurface> for BsplineSurface<Point3> {
    type Error = StepConvertingError;
    fn try_from(value: &NonRationalBsplineSurface) -> Result<Self, Self::Error> {
        use NonRationalBsplineSurface::{
            BezierSurface, BsplineSurfaceWithKnots, QuasiUniformSurface, UniformSurface,
        };
        match value {
            BsplineSurfaceWithKnots(x) => x.as_ref().try_into(),
            UniformSurface(x) => x.as_ref().try_into(),
            QuasiUniformSurface(x) => x.as_ref().try_into(),
            BezierSurface(x) => x.as_ref().try_into(),
        }
    }
}

/// `rational_b_spline_surface` as complex entity
///
/// This struct is an ad hoc implementation that differs from the definition by EXPRESS:
/// in AP042, rationalized curves are defined as complex entities,
/// but here the surfaces before rationalization are held as internal variables.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = rational_b_spline_surface)]
#[holder(generate_deserialize)]
pub struct RationalBsplineSurface {
    #[holder(use_place_holder)]
    non_rational_b_spline_surface: NonRationalBsplineSurface,
    weights_data: Vec<Vec<f64>>,
}

impl TryFrom<&RationalBsplineSurface> for NurbsSurface<Vector4> {
    type Error = StepConvertingError;
    fn try_from(
        RationalBsplineSurface {
            non_rational_b_spline_surface,
            weights_data,
        }: &RationalBsplineSurface,
    ) -> Result<Self, Self::Error> {
        let surface: BsplineSurface<Point3> = non_rational_b_spline_surface.try_into()?;
        Ok(Self::try_from_bspline_and_weights(
            surface,
            weights_data.clone(),
        )?)
    }
}

/// `swept_surface`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(generate_deserialize)]
pub enum SweptSurfaceAny {
    #[holder(use_place_holder)]
    SurfaceOfLinearExtrusion(Box<SurfaceOfLinearExtrusion>),
    #[holder(use_place_holder)]
    SurfaceOfRevolution(Box<SurfaceOfRevolution>),
}

impl TryFrom<&SweptSurfaceAny> for SweepSurface {
    type Error = StepConvertingError;
    fn try_from(value: &SweptSurfaceAny) -> Result<Self, Self::Error> {
        use SweptSurfaceAny::{SurfaceOfLinearExtrusion, SurfaceOfRevolution};
        Ok(match value {
            SurfaceOfLinearExtrusion(x) => Self::ExtrusionSurface(x.as_ref().try_into()?),
            SurfaceOfRevolution(x) => Self::RevolutionSurface(x.as_ref().try_into()?),
        })
    }
}

/// `surface_of_linear_extrusion`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = surface_of_linear_extrusion)]
#[holder(generate_deserialize)]
pub struct SurfaceOfLinearExtrusion {
    label: String,
    #[holder(use_place_holder)]
    swept_curve: CurveAny,
    #[holder(use_place_holder)]
    extrusion_axis: Vector,
}

impl TryFrom<&SurfaceOfLinearExtrusion> for StepExtrusionSurface {
    type Error = StepConvertingError;
    fn try_from(sr: &SurfaceOfLinearExtrusion) -> Result<Self, Self::Error> {
        let curve = Curve3D::try_from(&sr.swept_curve)?;
        let vector = Vector3::from(&sr.extrusion_axis);
        Ok(Self::by_extrusion(curve, vector))
    }
}

/// `surface_of_revolution`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = surface_of_revolution)]
#[holder(generate_deserialize)]
pub struct SurfaceOfRevolution {
    label: String,
    #[holder(use_place_holder)]
    swept_curve: CurveAny,
    #[holder(use_place_holder)]
    axis_position: Axis1Placement,
}

impl TryFrom<&SurfaceOfRevolution> for StepRevolutionSurface {
    type Error = StepConvertingError;
    fn try_from(sr: &SurfaceOfRevolution) -> Result<Self, Self::Error> {
        let curve = Curve3D::try_from(&sr.swept_curve)?;
        let origin = Point3::from(&sr.axis_position.location);
        let axis = sr.axis_position.direction().normalize();
        let mut rev = Self::new(RevolutionSurface::by_revolution(curve, origin, axis));
        rev.invert();
        Ok(rev)
    }
}
