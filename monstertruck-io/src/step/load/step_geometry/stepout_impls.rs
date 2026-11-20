use super::{
    ConicalSurface, Curve3D, ElementarySurface, InnerSpace, Matrix4, Point3, Processor, Sphere,
    SurfaceCurve3D, SurfaceCurveAssociatedGeometry, SurfaceCurveKind, SurfaceCurveRepresentation,
    Transform, Transformed, Vector3, save,
};
use crate::step::load::{GeometricCurveSet, GeometricSetSelect};
use crate::step::save::{
    FloatDisplay, IndexSliceDisplay, StepDisplay, StepLength, VectorAsDirection,
};

impl save::StepFormat for SurfaceCurveAssociatedGeometry {
    fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ParameterCurve(curve) => curve.fmt(idx, f),
            Self::Surface(surface) => surface.fmt(idx, f),
        }
    }
}

impl save::StepLength for SurfaceCurveAssociatedGeometry {
    fn step_length(&self) -> usize {
        match self {
            Self::ParameterCurve(curve) => curve.step_length(),
            Self::Surface(surface) => surface.step_length(),
        }
    }
}

impl save::StepFormat for SurfaceCurve3D {
    fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let leader_idx = idx + 1;
        let (associated_indices, _) = self.associated_geometry().iter().fold(
            (
                Vec::<usize>::with_capacity(self.associated_geometry().len()),
                leader_idx + self.leader().step_length(),
            ),
            |(mut indices, cursor), entry| {
                indices.push(cursor);
                (indices, cursor + entry.step_length())
            },
        );
        let entity = match self.kind() {
            SurfaceCurveKind::Surface => "SURFACE_CURVE",
            SurfaceCurveKind::Seam => "SEAM_CURVE",
            SurfaceCurveKind::Intersection => "INTERSECTION_CURVE",
        };
        let master_representation = match self.master_representation() {
            SurfaceCurveRepresentation::Curve3D => ".CURVE_3D.",
            SurfaceCurveRepresentation::ParameterCurve0 => ".PCURVE_S1.",
            SurfaceCurveRepresentation::ParameterCurve1 => ".PCURVE_S2.",
        };
        f.write_fmt(format_args!(
            "#{idx} = {entity}('', #{leader_idx}, {associated_geometry}, {master_representation});\n",
            associated_geometry = IndexSliceDisplay(associated_indices.iter().copied()),
        ))?;
        self.leader().fmt(leader_idx, f)?;
        self.associated_geometry()
            .iter()
            .zip(associated_indices)
            .try_for_each(|(entry, entry_idx)| entry.fmt(entry_idx, f))
    }
}

impl save::StepLength for SurfaceCurve3D {
    fn step_length(&self) -> usize {
        1 + self.leader().step_length()
            + self
                .associated_geometry()
                .iter()
                .map(save::StepLength::step_length)
                .sum::<usize>()
    }
}

impl save::StepCurve for SurfaceCurve3D {
    fn same_sense(&self) -> bool { self.leader().same_sense() }
}

impl save::ConstStepLength for Processor<Sphere, Matrix4> {
    const LENGTH: usize = Processor::<monstertruck_geometry::prelude::Sphere, Matrix4>::LENGTH;
}
impl save::StepLength for Processor<Sphere, Matrix4> {
    fn step_length(&self) -> usize { <Self as save::ConstStepLength>::LENGTH }
}
impl save::StepFormat for Processor<Sphere, Matrix4> {
    fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Processor::new(self.entity().0)
            .transformed(*self.transform())
            .fmt(idx, f)
    }
}

/// An exact STEP circular cone has a right-circular, not sheared or elliptic,
/// cross-section. Its radius and slope must be measured *after* applying
/// the processor's affine placement, including a length-unit transform.
struct StepConeFrame {
    location: Point3,
    axis: Vector3,
    reference: Vector3,
    radius: f64,
    semi_angle: f64,
}

fn source_cone_step_frame(processor: &ConicalSurface) -> Option<StepConeFrame> {
    let revo = processor.entity();
    let transform = processor.transform();
    let line = revo.entity_curve();
    let raw_axis = revo.axis();
    let axis_len = raw_axis.magnitude();
    if !axis_len.is_finite() || axis_len <= 0.0 {
        return None;
    }
    let local_axis = raw_axis / axis_len;
    let offset = line.0 - revo.origin();
    let local_foot = revo.origin() + local_axis * offset.dot(local_axis);
    let radial = line.0 - local_foot;
    let local_radius = radial.magnitude();
    if !local_radius.is_finite() || local_radius <= 0.0 {
        return None;
    }
    let local_reference = radial / local_radius;
    let local_tangent = local_axis.cross(local_reference);
    let derivative = line.1 - line.0;
    let radial_slope = derivative.dot(local_reference);
    let axial_slope = derivative.dot(local_axis);
    let tangential_slope = derivative.dot(local_tangent);
    let eps = 1.0e-11;
    if !axial_slope.is_finite()
        || axial_slope.abs() <= eps * derivative.magnitude()
        || tangential_slope.abs() > eps * derivative.magnitude()
    {
        // An oblique generator represents a different revolution surface.
        return None;
    }

    let world_reference = transform.transform_vector(local_reference);
    let world_tangent = transform.transform_vector(local_tangent);
    let world_axis = transform.transform_vector(local_axis);
    let rscale = world_reference.magnitude();
    let tscale = world_tangent.magnitude();
    let zscale = world_axis.magnitude();
    if [rscale, tscale, zscale]
        .into_iter()
        .any(|value| !value.is_finite() || value <= 0.0)
    {
        return None;
    }
    let orthogonal = |a: Vector3, na: f64, b: Vector3, nb: f64| (a.dot(b) / (na * nb)).abs() <= eps;
    if (rscale - tscale).abs() > eps * rscale.max(tscale)
        || !orthogonal(world_reference, rscale, world_tangent, tscale)
        || !orthogonal(world_reference, rscale, world_axis, zscale)
        || !orthogonal(world_tangent, tscale, world_axis, zscale)
    {
        // STEP's CONICAL_SURFACE cannot encode an elliptical or skew cone.
        return None;
    }

    let location = transform.transform_point(local_foot);
    let radius = local_radius * rscale;
    let semi_angle = (radial_slope * rscale / (axial_slope * zscale)).atan();
    if !radius.is_finite() || radius <= 0.0 || !semi_angle.is_finite() {
        return None;
    }
    Some(StepConeFrame {
        location,
        axis: world_axis,
        reference: world_reference,
        radius,
        semi_angle,
    })
}

impl save::StepFormat for ElementarySurface {
    fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Plane(x) => x.fmt(idx, f),
            Self::Sphere(x) => x.fmt(idx, f),
            Self::ToroidalSurface(x) => x.fmt(idx, f),
            Self::CylindricalSurface(processor) => {
                let position_idx = idx + 1;
                let location_idx = idx + 2;
                let axis_idx = idx + 3;
                let ref_direction_idx = idx + 4;

                let revo = processor.entity();
                let trans = processor.transform();
                let o = trans.transform_point(revo.origin());
                let p = trans.transform_point(revo.entity_curve().0);
                let axis = trans.transform_vector(revo.axis());

                let location = StepDisplay::new(o, location_idx);
                let direction_axis = VectorAsDirection(axis);
                let axis = StepDisplay::new(direction_axis, axis_idx);
                let raw_ref_direction = VectorAsDirection((p - o).normalize());
                let ref_direction = StepDisplay::new(raw_ref_direction, ref_direction_idx);
                let radius = (p - o).magnitude();

                f.write_fmt(format_args!(
                    "#{idx} = CYLINDRICAL_SURFACE('', #{position_idx}, {radius});
#{position_idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});
{location}{axis}{ref_direction}"
                ))
            }
            Self::ConicalSurface(processor) => {
                let frame = source_cone_step_frame(processor).ok_or(std::fmt::Error)?;
                let radius = FloatDisplay(frame.radius);
                let semi_angle = FloatDisplay(frame.semi_angle);

                let position_idx = idx + 1;
                let location_idx = idx + 2;
                let axis_idx = idx + 3;
                let ref_direction_idx = idx + 4;

                let location = StepDisplay::new(frame.location, location_idx);
                // STEP DIRECTION values are direction ratios, not
                // transformed basis vectors: never serialize their scale.
                let axis = StepDisplay::new(VectorAsDirection(frame.axis.normalize()), axis_idx);
                let ref_direction = StepDisplay::new(
                    VectorAsDirection(frame.reference.normalize()),
                    ref_direction_idx,
                );

                f.write_fmt(format_args!(
                    "#{idx} = CONICAL_SURFACE('', #{position_idx}, {radius}, {semi_angle});
#{position_idx} = AXIS2_PLACEMENT_3D('', #{location_idx}, #{axis_idx}, #{ref_direction_idx});
{location}{axis}{ref_direction}"
                ))
            }
        }
    }
}

/// Step length of a single [`GeometricSetSelect`] element.
fn set_select_step_length(elem: &GeometricSetSelect) -> usize {
    match elem {
        GeometricSetSelect::Curve(c) => {
            Curve3D::try_from(c.as_ref()).map_or(0, |c3d| c3d.step_length())
        }
        // Point3 has ConstStepLength = 1.
        GeometricSetSelect::Point(_) => 1,
    }
}

impl save::StepLength for GeometricCurveSet {
    fn step_length(&self) -> usize {
        // 1 for the GEOMETRIC_CURVE_SET entity + sum of element lengths.
        1 + self
            .elements
            .iter()
            .map(set_select_step_length)
            .sum::<usize>()
    }
}

impl save::StepFormat for GeometricCurveSet {
    fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Collect element start indices.
        let mut cursor = idx + 1;
        let element_indices: Vec<usize> = self
            .elements
            .iter()
            .filter_map(|e| {
                let len = set_select_step_length(e);
                if len == 0 {
                    return None;
                }
                let this = cursor;
                cursor += len;
                Some(this)
            })
            .collect();

        // Write the GEOMETRIC_CURVE_SET entity.
        f.write_fmt(format_args!(
            "#{idx} = GEOMETRIC_CURVE_SET('{}', {});\n",
            self.label,
            IndexSliceDisplay(element_indices.into_iter()),
        ))?;

        // Write each element.
        let mut cursor = idx + 1;
        for elem in &self.elements {
            match elem {
                GeometricSetSelect::Curve(c) => {
                    if let Ok(c3d) = Curve3D::try_from(c.as_ref()) {
                        save::StepFormat::fmt(&c3d, cursor, f)?;
                        cursor += c3d.step_length();
                    }
                }
                GeometricSetSelect::Point(p) => {
                    let pt = Point3::from(p.as_ref());
                    save::StepFormat::fmt(&pt, cursor, f)?;
                    cursor += 1;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod cone_serialization_tests {
    use super::*;
    use monstertruck_geometry::prelude::{Line, RevolutionSurface, Vector3};

    fn cone(p: Point3, direction: Vector3, transform: Matrix4) -> ConicalSurface {
        let local = RevolutionSurface::by_revolution(
            Line(p, p + direction),
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
        );
        Processor::with_transform(local, transform)
    }

    fn serialized_direction_lengths(fragment: &str) -> Vec<f64> {
        use step_p21::ast::{EntityInstance, Parameter};
        let document = format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('cone directions'),'2;1');\n\
FILE_NAME('cone','',(''),(''),'','','');\n\
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{fragment}\
ENDSEC;\nEND-ISO-10303-21;\n"
        );
        let parsed = step_p21::parser::parse(&document).expect("serialized cone is valid STEP");
        parsed.data[0]
            .entities
            .iter()
            .filter_map(|entity| {
                let EntityInstance::Simple { record, .. } = entity else {
                    return None;
                };
                if record.name != "DIRECTION" {
                    return None;
                }
                let Parameter::List(params) = &record.parameter else {
                    panic!("DIRECTION has no STEP parameter list")
                };
                let Some(Parameter::List(components)) = params.get(1) else {
                    panic!("DIRECTION is missing its coordinate ratios")
                };
                assert_eq!(components.len(), 3);
                let norm2 = components
                    .iter()
                    .map(|component| match component {
                        Parameter::Real(x) => x * x,
                        Parameter::Integer(x) => (*x as f64).powi(2),
                        _ => panic!("DIRECTION has a non-numeric ratio"),
                    })
                    .sum::<f64>();
                Some(norm2.sqrt())
            })
            .collect()
    }

    #[test]
    fn cone_radius_and_slope_follow_source_coordinate_scale() {
        let original = cone(
            Point3::new(0.3302, 0.0, 0.0),
            Vector3::new(0.6, 0.0, 1.0),
            Matrix4::from_scale(1.0),
        );
        let converted = cone(
            Point3::new(0.3302, 0.0, 0.0),
            Vector3::new(0.6, 0.0, 1.0),
            Matrix4::from_scale(0.1),
        );
        let orig = source_cone_step_frame(&original).expect("circular source cone");
        let scaled = source_cone_step_frame(&converted).expect("scaled circular source cone");
        assert!((orig.radius - 0.3302).abs() < 1.0e-14);
        assert!((scaled.radius - 0.03302).abs() < 1.0e-14);
        assert!((scaled.semi_angle - orig.semi_angle).abs() < 1.0e-14);
        let encoded = StepDisplay::new(ElementarySurface::ConicalSurface(converted), 1).to_string();
        assert!(encoded.contains("CONICAL_SURFACE"));
        assert!(encoded.contains("0.03302,"));
        assert!(!encoded.contains("0.3302,"));
        let lengths = serialized_direction_lengths(&encoded);
        assert_eq!(lengths.len(), 2);
        assert!(lengths.iter().all(|&length| (length - 1.0).abs() < 1.0e-12));
    }

    #[test]
    fn circular_cone_allows_distinct_axial_scale_and_general_radial_generator() {
        let transformed = cone(
            Point3::new(0.0, 0.3302, 0.4),
            Vector3::new(0.0, 0.6, 1.0),
            Matrix4::from_nonuniform_scale(0.1, 0.1, 0.2),
        );
        let frame = source_cone_step_frame(&transformed).expect("circular cone with axial stretch");
        assert!((frame.radius - 0.03302).abs() < 1.0e-14);
        assert!((frame.semi_angle - 0.3_f64.atan()).abs() < 1.0e-14);
        assert!((frame.location.z - 0.08).abs() < 1.0e-14);
    }

    #[test]
    fn source_export_rejects_elliptic_radial_scaling_and_oblique_generator() {
        let uneven = cone(
            Point3::new(0.3302, 0.0, 0.0),
            Vector3::new(0.6, 0.0, 1.0),
            Matrix4::from_nonuniform_scale(0.1, 0.2, 0.1),
        );
        assert!(source_cone_step_frame(&uneven).is_none());

        let skew = cone(
            Point3::new(0.3302, 0.0, 0.0),
            Vector3::new(0.6, 0.1, 1.0),
            Matrix4::from_scale(0.1),
        );
        assert!(source_cone_step_frame(&skew).is_none());
    }
}
