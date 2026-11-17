//! Exact complex-entity layouts for rational B-spline curve and surface imports.
//!
//! Each STEP subtype has seven sorted component records, but the position of
//! the weights, base B-spline parameters and label differs. Declare the mapping
//! once rather than repeating similar positional parsing in four dispatchers.
use super::{
    NonRationalBsplineCurveHolder as Curve, NonRationalBsplineSurfaceHolder as Surface,
    RationalBsplineCurveHolder, RationalBsplineSurfaceHolder, Table,
};
use serde::Deserialize;
use step_p21::{
    ast::{Parameter, Record},
    error::{Error, Result},
    tables::PlaceHolder,
};

#[derive(Clone, Copy)]
enum Geometry {
    Curve,
    Surface,
}

#[derive(Clone, Copy)]
enum Basis {
    Knots,
    Bezier,
    QuasiUniform,
    Uniform,
}

// STEP complex components are sorted alphabetically. The six shared
// components plus one of four supported subtype records define each
// recognized curve/surface layout exactly, without positional constants.
const CURVE_CORE: [&str; 6] = [
    "BOUNDED_CURVE",
    "B_SPLINE_CURVE",
    "CURVE",
    "GEOMETRIC_REPRESENTATION_ITEM",
    "RATIONAL_B_SPLINE_CURVE",
    "REPRESENTATION_ITEM",
];

const SURFACE_CORE: [&str; 6] = [
    "BOUNDED_SURFACE",
    "B_SPLINE_SURFACE",
    "GEOMETRIC_REPRESENTATION_ITEM",
    "RATIONAL_B_SPLINE_SURFACE",
    "REPRESENTATION_ITEM",
    "SURFACE",
];

fn identify_layout(records: &[Record]) -> Option<(Geometry, Basis)> {
    if records.len() != 7 {
        return None;
    }
    let (geometry, required) = if records
        .iter()
        .any(|record| record.name == "RATIONAL_B_SPLINE_CURVE")
    {
        (Geometry::Curve, &CURVE_CORE)
    } else if records
        .iter()
        .any(|record| record.name == "RATIONAL_B_SPLINE_SURFACE")
    {
        (Geometry::Surface, &SURFACE_CORE)
    } else {
        return None;
    };
    let mut extra = records
        .iter()
        .filter(|record| !required.contains(&record.name.as_str()));
    let subtype = extra.next()?.name.as_str();
    if extra.next().is_some() {
        return None;
    }
    let basis = match (geometry, subtype) {
        (Geometry::Curve, "B_SPLINE_CURVE_WITH_KNOTS")
        | (Geometry::Surface, "B_SPLINE_SURFACE_WITH_KNOTS") => Basis::Knots,
        (Geometry::Curve, "BEZIER_CURVE") | (Geometry::Surface, "BEZIER_SURFACE") => Basis::Bezier,
        (Geometry::Curve, "QUASI_UNIFORM_CURVE") | (Geometry::Surface, "QUASI_UNIFORM_SURFACE") => {
            Basis::QuasiUniform
        }
        (Geometry::Curve, "UNIFORM_CURVE") | (Geometry::Surface, "UNIFORM_SURFACE") => {
            Basis::Uniform
        }
        _ => return None,
    };
    // Preserve the original strict alphabetical component ordering. This
    // also rejects missing/duplicated base components and unknown subtypes.
    let mut required = required.iter().copied().peekable();
    let mut extra_used = false;
    for record in records {
        let expected = if !extra_used && required.peek().is_none_or(|name| subtype < *name) {
            extra_used = true;
            subtype
        } else {
            required.next()?
        };
        if record.name != expected {
            return None;
        }
    }
    if !extra_used || required.next().is_some() {
        return None;
    }
    Some((geometry, basis))
}

fn list_component<'a>(records: &'a [Record], name: &str) -> Option<&'a [Parameter]> {
    let record = records.iter().find(|record| record.name == name)?;
    match &record.parameter {
        Parameter::List(values) => Some(values),
        _ => None,
    }
}

struct Components<'a> {
    geometry: Geometry,
    basis: Basis,
    params: Parameter,
    weights: &'a Parameter,
}

fn unpack_components(records: &[Record]) -> Result<Option<Components<'_>>> {
    let Some((geometry, basis)) = identify_layout(records) else {
        return Ok(None);
    };
    let (base_name, weight_name, knot_name) = match geometry {
        Geometry::Curve => (
            "B_SPLINE_CURVE",
            "RATIONAL_B_SPLINE_CURVE",
            "B_SPLINE_CURVE_WITH_KNOTS",
        ),
        Geometry::Surface => (
            "B_SPLINE_SURFACE",
            "RATIONAL_B_SPLINE_SURFACE",
            "B_SPLINE_SURFACE_WITH_KNOTS",
        ),
    };
    let (Some(base), Some(label), Some(weights)) = (
        list_component(records, base_name),
        list_component(records, "REPRESENTATION_ITEM"),
        list_component(records, weight_name),
    ) else {
        return Ok(None);
    };
    let mut params = Vec::with_capacity(base.len() + label.len() + 5);
    match basis {
        Basis::Knots | Basis::Bezier => params.extend(label.iter().cloned()),
        Basis::QuasiUniform | Basis::Uniform => params.push(
            label
                .first()
                .ok_or_else(|| {
                    Error::DeserializeFailed("rational STEP entity has no label".to_owned())
                })?
                .clone(),
        ),
    }
    params.extend(base.iter().cloned());
    if let Basis::Knots = basis {
        let Some(knots) = list_component(records, knot_name) else {
            return Ok(None);
        };
        params.extend(knots.iter().cloned());
    }
    let weight = weights.first().ok_or_else(|| {
        Error::DeserializeFailed("rational STEP entity has no weight data".to_owned())
    })?;
    Ok(Some(Components {
        geometry,
        basis,
        params: Parameter::List(params),
        weights: weight,
    }))
}

/// Curve and surface non-rational enums have matching basis families
/// but different payload types; share the dispatch while retaining typed constructors.
macro_rules! decode_basis {
    ($basis:expr, $params:expr, $owner:ident, $knots:ident, $bezier:ident, $quasi:ident, $uniform:ident) => {
        match $basis {
            Basis::Knots => $owner::$knots(Deserialize::deserialize(&$params)?),
            Basis::Bezier => $owner::$bezier(Deserialize::deserialize(&$params)?),
            Basis::QuasiUniform => $owner::$quasi(Deserialize::deserialize(&$params)?),
            Basis::Uniform => $owner::$uniform(Deserialize::deserialize(&$params)?),
        }
    };
}

impl Table {
    pub(super) fn push_rational_complex(&mut self, id: u64, records: &[Record]) -> Result<bool> {
        let Some(Components {
            geometry,
            basis,
            params,
            weights,
        }) = unpack_components(records)?
        else {
            return Ok(false);
        };

        match geometry {
            Geometry::Curve => {
                let non_rational = decode_basis!(
                    basis,
                    params,
                    Curve,
                    BsplineCurveWithKnots,
                    BezierCurve,
                    QuasiUniformCurve,
                    UniformCurve
                );
                self.rational_b_spline_curve.insert(
                    id,
                    RationalBsplineCurveHolder {
                        non_rational_b_spline_curve: PlaceHolder::Owned(non_rational),
                        weights_data: Deserialize::deserialize(weights)?,
                    },
                );
            }
            Geometry::Surface => {
                let non_rational = decode_basis!(
                    basis,
                    params,
                    Surface,
                    BsplineSurfaceWithKnots,
                    BezierSurface,
                    QuasiUniformSurface,
                    UniformSurface
                );
                self.rational_b_spline_surface.insert(
                    id,
                    RationalBsplineSurfaceHolder {
                        non_rational_b_spline_surface: PlaceHolder::Owned(non_rational),
                        weights_data: Deserialize::deserialize(weights)?,
                    },
                );
            }
        }
        Ok(true)
    }
}
