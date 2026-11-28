//! Source-declared STEP length uncertainty and the effective recovery budget.
//! Do not use a candidate export's metadata to relax its source comparison.
use crate::step_entities::{parameter_number, representation_items_and_context};
use crate::step_graph::{entity_ref_value, simple_record};
use monstertruck_core::tolerance::{NATIVE_STEP_RECOVERY_MM, native_step_recovery_mm};
use std::collections::HashMap;
use step_p21::ast::{EntityInstance, Parameter, Record};

/// Effective native recovery budgets for representations declaring uncertainty.
/// A missing entry uses the project default; author uncertainty is a lower bound.
pub fn source_tolerance_by_representation_item(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> HashMap<u64, f64> {
    declared_source_uncertainty_by_representation_item(entities, index)
        .into_iter()
        .map(|(item, uncertainty)| (item, native_step_recovery_mm(Some(uncertainty))))
        .collect()
}

/// Effective tolerance for one source body (never below its declared precision).
pub fn source_recovery_tolerance_mm(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> f64 {
    source_tolerance_by_representation_item(entities, index)
        .get(&solid_id)
        .copied()
        .unwrap_or(NATIVE_STEP_RECOVERY_MM)
}

/// Literal author-declared length uncertainty (including unit conversion),
/// without any fallback to the kernel's own recovery geometry tolerance.
pub fn declared_source_uncertainty_by_representation_item(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> HashMap<u64, f64> {
    let mut context_cache = HashMap::<u64, Option<f64>>::new();
    let mut out = HashMap::<u64, f64>::new();

    for entity in entities {
        let Some(record) = simple_record(entity) else {
            continue;
        };
        if !(record.name == "REPRESENTATION" || record.name.ends_with("_REPRESENTATION")) {
            continue;
        }
        let Some((items, context_id)) = representation_items_and_context(entity) else {
            continue;
        };
        let tolerance = *context_cache.entry(context_id).or_insert_with(|| {
            representation_context_length_uncertainty_mm(context_id, entities, index)
        });
        let Some(tolerance) = tolerance else {
            continue;
        };

        for item in items {
            out.entry(item)
                .and_modify(|existing| *existing = existing.max(tolerance))
                .or_insert(tolerance);
        }
    }

    out
}

fn representation_context_length_uncertainty_mm(
    context_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<f64> {
    let entity = entities.get(*index.get(&context_id)?)?;
    let record = entity_record_named(entity, "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT")?;
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    let [Parameter::List(uncertainties)] = params.as_slice() else {
        return None;
    };

    uncertainties
        .iter()
        .filter_map(entity_ref_value)
        .filter_map(|id| uncertainty_measure_mm(id, entities, index))
        .filter(|value| value.is_finite() && *value > 0.0)
        .reduce(f64::max)
}

fn uncertainty_measure_mm(
    uncertainty_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<f64> {
    let entity = entities.get(*index.get(&uncertainty_id)?)?;
    let record = entity_record_named(entity, "UNCERTAINTY_MEASURE_WITH_UNIT")?;
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    let [
        Parameter::Typed { keyword, parameter },
        unit,
        Parameter::String(_),
        _,
    ] = params.as_slice()
    else {
        return None;
    };
    if keyword.as_ref() != "LENGTH_MEASURE" {
        return None;
    }
    let value = parameter_number(parameter)?;
    let unit_id = entity_ref_value(unit)?;
    let section = crate::step_document::StepSectionView::new(entities, index);
    let scale_mm = crate::units::length_unit_scale_mm(unit_id, &section)?;
    let result = value * scale_mm;
    result.is_finite().then_some(result)
}

fn entity_record_named<'a>(entity: &'a EntityInstance, name: &str) -> Option<&'a Record> {
    match entity {
        EntityInstance::Simple { record, .. } => (record.name == name).then_some(record),
        EntityInstance::Complex { subsuper, .. } => {
            subsuper.0.iter().find(|record| record.name == name)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_graph::{build_index, entity_id};
    #[test]
    fn reads_representation_length_uncertainty_without_an_upper_cap() -> anyhow::Result<()> {
        fn tolerance_from(text: &str) -> anyhow::Result<(Option<f64>, f64)> {
            let exchange = step_p21::parser::parse(text)?;
            let entities = &exchange
                .data
                .first()
                .ok_or_else(|| anyhow::anyhow!("parsed exchange has no DATA section"))?
                .entities;
            let index = build_index(entities);
            let solid_id = entities
                .iter()
                .find(|entity| {
                    simple_record(entity).is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
                })
                .map(entity_id)
                .ok_or_else(|| anyhow::anyhow!("fixture has no MANIFOLD_SOLID_BREP"))?;
            Ok((
                declared_source_uncertainty_by_representation_item(entities, &index)
                    .get(&solid_id)
                    .copied(),
                source_recovery_tolerance_mm(solid_id, entities, &index),
            ))
        }
        let source = std::str::from_utf8(include_bytes!(
            "../validation/fixtures/native_conical_frustum.step"
        ))?;
        assert_eq!(tolerance_from(source)?, (Some(1.0e-7), 0.01));
        let widened = source.replacen("LENGTH_MEASURE(1.E-07)", "LENGTH_MEASURE(5.E-02)", 1);
        assert_eq!(tolerance_from(&widened)?, (Some(0.05), 0.05));
        let centimetres = widened.replacen(
            "LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.)",
            "LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.CENTI.,.METRE.)",
            1,
        );
        assert_eq!(tolerance_from(&centimetres)?, (Some(0.5), 0.5));
        let duplicate = widened.replacen(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117,#117))",
            1,
        );
        assert_eq!(tolerance_from(&duplicate)?, (Some(0.05), 0.05));
        let absent = source.replace(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(())",
        );
        assert_eq!(tolerance_from(&absent)?, (None, 0.01));
        Ok(())
    }
}
