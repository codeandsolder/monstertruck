//! Regression coverage for the eight recognized rational complex STEP layouts.
//! These tests exercise the real typed holder deserializers, not only record-name matching.
use super::*;
use std::str::FromStr;

const CURVE_BSPLINE: &str = "B_SPLINE_CURVE(2,(#1,#1,#1),.UNSPECIFIED.,.U.,.U.)";
const SURFACE_BSPLINE: &str =
    "B_SPLINE_SURFACE(2,2,((#1,#1,#1),(#1,#1,#1),(#1,#1,#1)),.UNSPECIFIED.,.U.,.U.,.U.)";
const CURVE_WEIGHTS: &str = "RATIONAL_B_SPLINE_CURVE((1.0,2.0,3.0))";
const SURFACE_WEIGHTS: &str =
    "RATIONAL_B_SPLINE_SURFACE(((1.0,2.0,3.0),(4.0,5.0,6.0),(7.0,8.0,9.0)))";

fn rational_instance(surface: bool, basis: &str) -> EntityInstance {
    let (bspline, weights, bounded, base, base_curve, representation) = if surface {
        (
            SURFACE_BSPLINE,
            SURFACE_WEIGHTS,
            "BOUNDED_SURFACE()",
            "SURFACE()",
            "GEOMETRIC_REPRESENTATION_ITEM()",
            "REPRESENTATION_ITEM('surface')",
        )
    } else {
        (
            CURVE_BSPLINE,
            CURVE_WEIGHTS,
            "BOUNDED_CURVE()",
            "CURVE()",
            "GEOMETRIC_REPRESENTATION_ITEM()",
            "REPRESENTATION_ITEM('curve')",
        )
    };
    let extra = match (surface, basis) {
        (false, "knots") => "B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.0,1.0),.UNSPECIFIED.)",
        (true, "knots") => {
            "B_SPLINE_SURFACE_WITH_KNOTS((3,3),(3,3),(0.0,1.0),(0.0,1.0),.UNSPECIFIED.)"
        }
        (false, "bezier") => "BEZIER_CURVE()",
        (true, "bezier") => "BEZIER_SURFACE()",
        (false, "quasi") => "QUASI_UNIFORM_CURVE()",
        (true, "quasi") => "QUASI_UNIFORM_SURFACE()",
        (false, "uniform") => "UNIFORM_CURVE()",
        (true, "uniform") => "UNIFORM_SURFACE()",
        _ => unreachable!(),
    };
    let mut records = vec![
        bounded,
        bspline,
        extra,
        base,
        base_curve,
        weights,
        representation,
    ];
    if surface {
        records = vec![
            bounded,
            bspline,
            extra,
            base_curve,
            weights,
            representation,
            base,
        ];
    }
    records.sort_by_key(|record| record.split('(').next().expect("valid test record"));
    let source = format!("#101=({});", records.join(" "));
    EntityInstance::from_str(&source).expect("valid complex STEP record")
}

#[test]
fn all_eight_rational_complex_layouts_preserve_weights_and_kind() {
    for surface in [false, true] {
        for basis in ["knots", "bezier", "quasi", "uniform"] {
            let entity = rational_instance(surface, basis);
            let mut table = Table::default();
            table
                .push_instance(&entity)
                .expect("recognized rational layout");
            assert!(table.dummy.is_empty(), "{surface:?} {basis}");
            if surface {
                let holder = table
                    .rational_b_spline_surface
                    .get(&101)
                    .unwrap_or_else(|| panic!("missing {basis} surface"));
                assert!(format!("{holder:?}").contains("weights_data"));
                assert!(table.rational_b_spline_curve.is_empty());
            } else {
                let holder = table
                    .rational_b_spline_curve
                    .get(&101)
                    .unwrap_or_else(|| panic!("missing {basis} curve"));
                assert_eq!(holder.weights_data, [1.0, 2.0, 3.0]);
                assert!(table.rational_b_spline_surface.is_empty());
            }
        }
    }
}

#[test]
fn malformed_rational_weights_are_rejected_without_inserting_a_holder() {
    for surface in [false, true] {
        let mut entity = rational_instance(surface, "knots");
        let EntityInstance::Complex {
            subsuper: SubSuperRecord(records),
            ..
        } = &mut entity
        else {
            unreachable!()
        };
        let record = records
            .iter_mut()
            .find(|record| record.name.starts_with("RATIONAL_B_SPLINE_"))
            .expect("rational subtype exists");
        record.parameter = Parameter::List(Vec::new());

        let mut table = Table::default();
        assert!(table.push_instance(&entity).is_err());
        assert!(table.rational_b_spline_curve.is_empty());
        assert!(table.rational_b_spline_surface.is_empty());
    }
}

#[test]
fn missing_label_or_unknown_rational_layout_never_creates_a_typed_holder() {
    for surface in [false, true] {
        for basis in ["quasi", "uniform"] {
            let mut entity = rational_instance(surface, basis);
            let EntityInstance::Complex {
                subsuper: SubSuperRecord(records),
                ..
            } = &mut entity
            else {
                unreachable!()
            };
            records
                .iter_mut()
                .find(|record| record.name == "REPRESENTATION_ITEM")
                .expect("label record present")
                .parameter = Parameter::List(Vec::new());

            let mut table = Table::default();
            assert!(table.push_instance(&entity).is_err(), "{surface:?} {basis}");
            assert!(table.rational_b_spline_curve.is_empty());
            assert!(table.rational_b_spline_surface.is_empty());
        }

        let mut entity = rational_instance(surface, "knots");
        let EntityInstance::Complex {
            subsuper: SubSuperRecord(records),
            ..
        } = &mut entity
        else {
            unreachable!()
        };
        records
            .iter_mut()
            .find(|record| record.name.starts_with("RATIONAL_B_SPLINE_"))
            .expect("rational subtype present")
            .name
            .push_str("_UNKNOWN");

        let mut table = Table::default();
        table
            .push_instance(&entity)
            .expect("unrecognized layout falls back to dummy");
        assert!(table.rational_b_spline_curve.is_empty());
        assert!(table.rational_b_spline_surface.is_empty());
        assert!(table.dummy.contains_key(&101));
    }
}
