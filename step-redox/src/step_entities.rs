use crate::step_graph::{entity_ref_value, simple_record, simple_record_mut, visit_entity_refs};
use std::collections::{HashMap, HashSet};
use step_p21::ast::{EntityInstance, Name, Parameter, Record};

/// Resolve a STEP reference only when the index still names the requested entity.
/// An in-bounds but stale index must not silently substitute another entity.
fn indexed_entity<'a>(
    id: u64,
    entities: &'a [EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<&'a EntityInstance> {
    let entity = entities.get(*index.get(&id)?)?;
    (crate::step_graph::entity_id(entity) == id).then_some(entity)
}

/// Fetch the positional parameters of a simple STEP entity of the requested type.
/// Complex entities, malformed parameter lists, and stale IDs fail closed.
pub fn named_entity_parameters<'a>(
    id: u64,
    expected_name: &str,
    entities: &'a [EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<&'a [Parameter]> {
    let record = simple_record(indexed_entity(id, entities, index)?)?;
    if record.name != expected_name {
        return None;
    }
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    Some(params)
}

/// Read an exactly three-component STEP coordinate tuple of the required entity kind.
fn coordinate_components(
    id: u64,
    name: &str,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<[f64; 3]> {
    let params = named_entity_parameters(id, name, entities, index)?;
    let Parameter::List(components) = params.get(1)? else {
        return None;
    };
    let [x, y, z] = components.as_slice() else {
        return None;
    };
    Some([number(x)?, number(y)?, number(z)?])
}

pub fn cartesian_point(
    id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<[f64; 3]> {
    coordinate_components(id, "CARTESIAN_POINT", entities, index)
}

pub fn number(param: &Parameter) -> Option<f64> {
    match param {
        Parameter::Real(value) => Some(*value),
        Parameter::Integer(value) => crate::numeric::exact_i64_to_f64(*value),
        _ => None,
    }
}

pub fn direction_components(
    id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<[f64; 3]> {
    coordinate_components(id, "DIRECTION", entities, index)
}

pub fn vertex_point(
    vertex_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<[f64; 3]> {
    let record = simple_record(indexed_entity(vertex_id, entities, index)?)?;
    if record.name != "VERTEX_POINT" {
        return None;
    }
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    cartesian_point(entity_ref_value(params.get(1)?)?, entities, index)
}

pub fn enumeration_bool(parameter: &Parameter) -> Option<bool> {
    match parameter {
        Parameter::Enumeration(value) if value == "T" => Some(true),
        Parameter::Enumeration(value) if value == "F" => Some(false),
        _ => None,
    }
}

pub fn nth_entity_ref(parameter: &Parameter, idx: usize) -> Option<u64> {
    let Parameter::List(params) = parameter else {
        return None;
    };
    entity_ref_value(params.get(idx)?)
}

pub fn closure_from(
    root: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> HashSet<u64> {
    let mut seen = HashSet::new();
    let mut stack = vec![root];
    while let Some(id) = stack.pop() {
        if !seen.insert(id) {
            continue;
        }
        let Some(entity) = indexed_entity(id, entities, index) else {
            continue;
        };
        visit_entity_refs(entity, &mut |child| {
            if indexed_entity(child, entities, index).is_some() && !seen.contains(&child) {
                stack.push(child);
            }
        });
    }
    seen
}

pub fn representation_items_and_context(entity: &EntityInstance) -> Option<(Vec<u64>, u64)> {
    let record = simple_record(entity)?;
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    if params.len() != 3 {
        return None;
    }
    let Parameter::List(items) = &params[1] else {
        return None;
    };
    let item_ids: Option<Vec<u64>> = items.iter().map(entity_ref_value).collect();
    Some((item_ids?, entity_ref_value(&params[2])?))
}

pub fn entity_ref_list(parameter: &Parameter) -> Option<Vec<u64>> {
    let Parameter::List(items) = parameter else {
        return None;
    };
    items.iter().map(entity_ref_value).collect()
}

pub const fn integer_value(parameter: &Parameter) -> Option<i64> {
    match parameter {
        Parameter::Integer(value) => Some(*value),
        _ => None,
    }
}

pub fn integer_list(parameter: &Parameter) -> Option<Vec<i64>> {
    let Parameter::List(items) = parameter else {
        return None;
    };
    items.iter().map(integer_value).collect()
}

pub fn numeric_list(parameter: &Parameter) -> Option<Vec<f64>> {
    let Parameter::List(items) = parameter else {
        return None;
    };
    items.iter().map(number).collect()
}

pub fn finite_numeric_list(parameter: &Parameter) -> Option<Vec<f64>> {
    let Parameter::List(items) = parameter else {
        return None;
    };
    items
        .iter()
        .map(|item| number(item).filter(|value| value.is_finite()))
        .collect()
}

pub fn parameter_number(parameter: &Parameter) -> Option<f64> {
    match parameter {
        Parameter::Typed { parameter, .. } => parameter_number(parameter),
        _ => number(parameter),
    }
}

pub fn edge_vertices(
    edge: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<[u64; 2]> {
    let params = named_entity_parameters(edge, "EDGE_CURVE", entities, index)?;
    Some([
        entity_ref_value(params.get(1)?)?,
        entity_ref_value(params.get(2)?)?,
    ])
}

pub fn oriented_edge_element(record: &Record) -> Option<u64> {
    if record.name != "ORIENTED_EDGE" {
        return None;
    }
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    entity_ref_value(params.get(3)?)
}

pub fn oriented_edge_orientation(record: &Record) -> Option<bool> {
    let Parameter::List(params) = &record.parameter else {
        return None;
    };
    enumeration_bool(params.get(4)?)
}

pub fn styled_items_by_target(entities: &[EntityInstance]) -> HashMap<u64, Vec<(u64, Vec<u64>)>> {
    let mut out: HashMap<u64, Vec<(u64, Vec<u64>)>> = HashMap::new();
    for entity in entities {
        let id = crate::step_graph::entity_id(entity);
        let Some(record) = simple_record(entity) else {
            continue;
        };
        if record.name != "STYLED_ITEM" {
            continue;
        }
        let Parameter::List(params) = &record.parameter else {
            continue;
        };
        if params.len() != 3 {
            continue;
        }
        let Some(target) = entity_ref_value(&params[2]) else {
            continue;
        };
        let Some(assignments) = entity_ref_list(&params[1]) else {
            continue;
        };
        out.entry(target).or_default().push((id, assignments));
    }
    out
}

pub fn patch_presentation_lists(
    entities: &mut [EntityInstance],
    remove: &HashSet<u64>,
    add: &[u64],
) {
    for entity in entities {
        let Some(record) = simple_record_mut(entity) else {
            continue;
        };
        let target_index = match record.name.as_str() {
            "MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION" => 1,
            "PRESENTATION_LAYER_ASSIGNMENT" => 2,
            _ => continue,
        };
        let Parameter::List(params) = &mut record.parameter else {
            continue;
        };
        let Some(Parameter::List(items)) = params.get_mut(target_index) else {
            continue;
        };
        let had_removed = items
            .iter()
            .filter_map(entity_ref_value)
            .any(|id| remove.contains(&id));
        if !had_removed {
            continue;
        }
        items.retain(|item| entity_ref_value(item).is_none_or(|id| !remove.contains(&id)));
        items.extend(add.iter().copied().map(entity_ref));
    }
}

pub fn push_point(entities: &mut Vec<EntityInstance>, next_id: &mut u64, point: [f64; 3]) -> u64 {
    push_simple(
        entities,
        next_id,
        "CARTESIAN_POINT",
        vec![
            Parameter::String(String::new()),
            Parameter::List(point.into_iter().map(Parameter::Real).collect()),
        ],
    )
}

pub fn push_simple(
    entities: &mut Vec<EntityInstance>,
    next_id: &mut u64,
    name: &str,
    params: Vec<Parameter>,
) -> u64 {
    let id = *next_id;
    *next_id += 1;
    entities.push(EntityInstance::Simple {
        id,
        record: Record {
            name: name.to_string(),
            parameter: Parameter::List(params),
        },
    });
    id
}

pub const fn entity_ref(id: u64) -> Parameter { Parameter::Ref(Name::Entity(id)) }

#[cfg(test)]
mod coordinate_lookup_tests {
    use super::*;

    fn entity(id: u64, kind: &str, values: Vec<Parameter>) -> EntityInstance {
        EntityInstance::Simple {
            id,
            record: Record {
                name: kind.to_owned(),
                parameter: Parameter::List(vec![
                    Parameter::String(String::new()),
                    Parameter::List(values),
                ]),
            },
        }
    }

    #[test]
    fn shared_coordinate_accessors_preserve_point_and_direction_types() {
        let entities = vec![
            entity(
                1,
                "CARTESIAN_POINT",
                vec![
                    Parameter::Integer(3),
                    Parameter::Real(-2.5),
                    Parameter::Real(0.0),
                ],
            ),
            entity(
                2,
                "DIRECTION",
                vec![
                    Parameter::Real(0.0),
                    Parameter::Integer(0),
                    Parameter::Real(1.0),
                ],
            ),
        ];
        let index = HashMap::from([(1, 0), (2, 1)]);
        assert_eq!(
            cartesian_point(1, &entities, &index),
            Some([3.0, -2.5, 0.0])
        );
        assert_eq!(
            direction_components(2, &entities, &index),
            Some([0.0, 0.0, 1.0])
        );
        assert_eq!(cartesian_point(2, &entities, &index), None);
        assert_eq!(direction_components(1, &entities, &index), None);
    }

    #[test]
    fn shared_coordinate_accessors_reject_incomplete_or_stale_references() {
        let entities = vec![
            entity(1, "CARTESIAN_POINT", vec![Parameter::Real(1.0)]),
            entity(
                2,
                "DIRECTION",
                vec![
                    Parameter::Real(1.0),
                    Parameter::String("bad".to_owned()),
                    Parameter::Real(3.0),
                ],
            ),
        ];
        let index = HashMap::from([(1, 0), (2, 1), (3, 100)]);
        assert_eq!(cartesian_point(1, &entities, &index), None);
        assert_eq!(direction_components(2, &entities, &index), None);
        assert_eq!(cartesian_point(3, &entities, &index), None);
        assert_eq!(direction_components(3, &entities, &index), None);
        assert_eq!(cartesian_point(4, &entities, &index), None);
    }

    #[test]
    fn in_bounds_stale_indices_cannot_substitute_other_entities() {
        let entities = vec![
            entity(
                1,
                "CARTESIAN_POINT",
                vec![
                    Parameter::Real(1.0),
                    Parameter::Real(2.0),
                    Parameter::Real(3.0),
                ],
            ),
            entity(
                2,
                "CARTESIAN_POINT",
                vec![
                    Parameter::Real(4.0),
                    Parameter::Real(5.0),
                    Parameter::Real(6.0),
                ],
            ),
            EntityInstance::Simple {
                id: 10,
                record: Record {
                    name: "VERTEX_POINT".to_owned(),
                    parameter: Parameter::List(vec![
                        Parameter::String(String::new()),
                        entity_ref(1),
                    ]),
                },
            },
            EntityInstance::Simple {
                id: 20,
                record: Record {
                    name: "EDGE_CURVE".to_owned(),
                    parameter: Parameter::List(vec![
                        Parameter::String(String::new()),
                        entity_ref(10),
                        entity_ref(10),
                    ]),
                },
            },
        ];
        let mut index = HashMap::from([
            (1, 0),
            (2, 1),
            (10, 2),
            (20, 3),
            (99, 0),
            (110, 2),
            (120, 3),
        ]);
        assert_eq!(cartesian_point(99, &entities, &index), None);
        assert_eq!(vertex_point(110, &entities, &index), None);
        assert_eq!(edge_vertices(120, &entities, &index), None);
        assert_eq!(vertex_point(10, &entities, &index), Some([1.0, 2.0, 3.0]));
        assert_eq!(edge_vertices(20, &entities, &index), Some([10, 10]));
        assert_eq!(
            closure_from(20, &entities, &index),
            HashSet::from([1, 10, 20])
        );

        // A stale index that points at a different, but equally typed, point
        // must fail closed rather than quietly supplying the wrong coordinates.
        index.insert(1, 1);
        assert_eq!(cartesian_point(1, &entities, &index), None);
        assert_eq!(vertex_point(10, &entities, &index), None);
        assert_eq!(closure_from(20, &entities, &index), HashSet::from([10, 20]));
    }
}
