use crate::{
    shape_identity,
    step_document::{StepEntity, StepSection, StepValue},
    step_entities, step_graph, units,
};
use anyhow::{Context, Result, bail, ensure};
use monstertruck_io::step::load::ParsedStepSource;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::hash::Hasher;
use std::path::Path;
use step_p21::ast::{DataSection, EntityInstance, Exchange, Parameter, Record};

pub mod batch;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanReport {
    manifold_solids: usize,
    analyzed_solids: usize,
    skipped_solids: usize,
    fallback_identity_solids: usize,
    solids_without_shape_representation: usize,
    solids: Vec<SolidScan>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct SolidScan {
    data_section: usize,
    solid_id: u64,
    geometry_fingerprint: Option<String>,
    identity_kind: IdentityKind,
    unit_scale_mm: Option<f64>,
    center_source_units: Option<[f64; 3]>,
    canonical_quarter_turn: Option<u8>,
    closure_entities: usize,
    vertices: usize,
    edges: usize,
    oriented_edges: usize,
    faces: usize,
    has_shape_representation: bool,
}

impl SolidScan {
    const fn topology(&self) -> [usize; 4] {
        [self.vertices, self.edges, self.oriented_edges, self.faces]
    }
}

fn digest(bytes: &[u8]) -> String { hex::encode(Sha256::digest(bytes)) }

#[derive(Debug, Serialize, Deserialize)]
pub struct ExtractedBody {
    fingerprint: String,
    sha256: String,
    bytes: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum IdentityKind {
    Canonical,
    FileLocalFallback,
}

#[derive(Clone, Copy)]
struct Fnv64(u64);

impl Fnv64 {
    const fn seeded(seed: u64) -> Self { Self(seed) }
}

impl Hasher for Fnv64 {
    fn finish(&self) -> u64 { self.0 }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0100_0000_01b3);
        }
    }
}

fn feed_u64(hashers: &mut [Fnv64; 2], value: u64) {
    for hasher in hashers {
        hasher.write(&value.to_le_bytes());
    }
}

fn feed_i64(hashers: &mut [Fnv64; 2], value: i64) {
    for hasher in hashers {
        hasher.write(&value.to_le_bytes());
    }
}

fn feed_bytes(hashers: &mut [Fnv64; 2], value: &[u8]) {
    feed_u64(hashers, value.len() as u64);
    for hasher in hashers {
        hasher.write(value);
    }
}

fn fingerprint(value: &shape_identity::ShapeKey) -> String {
    let mut hashers = [
        Fnv64::seeded(0xcbf2_9ce4_8422_2325),
        Fnv64::seeded(0x8422_2325_cbf2_9ce4),
    ];
    for count in [
        value.vertices,
        value.edges,
        value.oriented_edges,
        value.faces,
    ] {
        feed_u64(&mut hashers, count as u64);
    }
    feed_u64(&mut hashers, value.points.len() as u64);
    for point in &value.points {
        for coordinate in point {
            feed_i64(&mut hashers, *coordinate);
        }
    }
    for geometry in [&value.edge_geometry, &value.face_geometry] {
        feed_u64(&mut hashers, geometry.len() as u64);
        for (name, scalars) in geometry {
            feed_bytes(&mut hashers, name.as_bytes());
            feed_u64(&mut hashers, scalars.len() as u64);
            for scalar in scalars {
                feed_i64(&mut hashers, *scalar);
            }
        }
    }
    feed_bytes(&mut hashers, value.topology.as_bytes());
    format!("{:016x}{:016x}", hashers[0].finish(), hashers[1].finish())
}

fn parse_exchange(bytes: &[u8]) -> Result<Exchange> { Ok(ParsedStepSource::parse(bytes)?.exchange) }

fn representation_ownership(entities: &[EntityInstance]) -> HashMap<u64, (u64, u64)> {
    let mut owners = HashMap::new();
    for entity in entities {
        let Some(record) = step_graph::simple_record(entity) else {
            continue;
        };
        if record.name != "ADVANCED_BREP_SHAPE_REPRESENTATION"
            && record.name != "SHAPE_REPRESENTATION"
        {
            continue;
        }
        let Some((items, context_id)) = step_entities::representation_items_and_context(entity)
        else {
            continue;
        };
        let representation_id = step_graph::entity_id(entity);
        for item in items {
            owners
                .entry(item)
                .or_insert((representation_id, context_id));
        }
    }
    owners
}

fn shape_definition_root_for_representation(
    representation_id: u64,
    entities: &[EntityInstance],
) -> Option<(u64, u64)> {
    entities.iter().find_map(|entity| {
        let record = step_graph::simple_record(entity)?;
        if record.name != "SHAPE_DEFINITION_REPRESENTATION" {
            return None;
        }
        let Parameter::List(params) = &record.parameter else {
            return None;
        };
        let [definition, representation] = params.as_slice() else {
            return None;
        };
        (step_graph::entity_ref_value(representation)? == representation_id).then_some((
            step_graph::entity_id(entity),
            step_graph::entity_ref_value(definition)?,
        ))
    })
}

/// Each isolated body is a new standalone product, regardless of how (or
/// whether) the source's representation was linked to a product definition.
///
/// Reusing source product roots can silently leave instance representations
/// orphaned: a product may point at an assembly representation containing a
/// separate instance, while the actual source solid is owned by a different
/// `ADVANCED_BREP_SHAPE_REPRESENTATION`. OCCT then transfers zero roots.
fn isolated_product_root(representation_id: u64, first_id: u64) -> Vec<EntityInstance> {
    let reference = step_entities::entity_ref;
    let make = |offset: u64, name: &str, parameters: Vec<Parameter>| EntityInstance::Simple {
        id: first_id + offset,
        record: Record {
            name: name.into(),
            parameter: Parameter::List(parameters),
        },
    };
    let string = |value: &str| Parameter::String(value.into());
    vec![
        make(0, "APPLICATION_CONTEXT", vec![string("mechanical design")]),
        make(
            1,
            "PRODUCT_CONTEXT",
            vec![
                string("part definition"),
                reference(first_id),
                string("mechanical"),
            ],
        ),
        make(
            2,
            "PRODUCT",
            vec![
                string("isolated body"),
                string("isolated body"),
                string(""),
                Parameter::List(vec![reference(first_id + 1)]),
            ],
        ),
        make(
            3,
            "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE",
            vec![
                string(""),
                string(""),
                reference(first_id + 2),
                Parameter::Enumeration("NOT_KNOWN".into()),
            ],
        ),
        make(
            4,
            "PRODUCT_DEFINITION_CONTEXT",
            vec![
                string("part definition"),
                reference(first_id),
                string("design"),
            ],
        ),
        make(
            5,
            "PRODUCT_DEFINITION",
            vec![
                string("isolated body"),
                string(""),
                reference(first_id + 3),
                reference(first_id + 4),
            ],
        ),
        make(
            6,
            "PRODUCT_DEFINITION_SHAPE",
            vec![string(""), string(""), reference(first_id + 5)],
        ),
        make(
            7,
            "SHAPE_DEFINITION_REPRESENTATION",
            vec![reference(first_id + 6), reference(representation_id)],
        ),
    ]
}

fn representation_length_scale_mm(context_id: u64, section: &StepSection<'_>) -> Option<f64> {
    let params = section
        .entity(context_id)?
        .record_named("GLOBAL_UNIT_ASSIGNED_CONTEXT")?
        .value()
        .list()?;
    if params.len() != 1 {
        return None;
    }
    let unit_values = params.get(0)?.list()?;
    let mut scales = unit_values
        .iter()
        .filter_map(StepValue::entity_ref)
        .filter_map(|unit_id| units::length_unit_scale_mm(unit_id, &section.view()));
    let scale = scales.next()?;
    scales.next().is_none().then_some(scale)
}

#[derive(Debug, Clone, Copy)]
struct BasicSolidStats {
    closure_entities: usize,
    vertices: usize,
    edges: usize,
    oriented_edges: usize,
    faces: usize,
}

fn solid_basic_stats(
    root: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> BasicSolidStats {
    let closure = shape_identity::semantic_solid_closure(root, entities, index)
        .unwrap_or_else(|| step_entities::closure_from(root, entities, index));
    let mut stats = BasicSolidStats {
        closure_entities: closure.len(),
        vertices: 0,
        edges: 0,
        oriented_edges: 0,
        faces: 0,
    };
    for id in closure {
        let Some(&entity_index) = index.get(&id) else {
            continue;
        };
        let Some(record) = step_graph::simple_record(&entities[entity_index]) else {
            continue;
        };
        match record.name.as_str() {
            "VERTEX_POINT" => stats.vertices += 1,
            "EDGE_CURVE" => stats.edges += 1,
            "ORIENTED_EDGE" => stats.oriented_edges += 1,
            "ADVANCED_FACE" => stats.faces += 1,
            _ => {}
        }
    }
    stats
}

/// Scan the STEP file using the production body-identity implementation.
///
/// # Errors
/// Returns an error if the file cannot be read, decoded, or parsed.
pub fn scan(path: &Path) -> Result<ScanReport> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    scan_bytes(&bytes)
}

fn scan_bytes(bytes: &[u8]) -> Result<ScanReport> {
    let exchange = parse_exchange(bytes)?;
    Ok(scan_exchange(&exchange))
}

fn scan_exchange(exchange: &Exchange) -> ScanReport {
    let mut solids = Vec::new();
    let mut manifold_solids = 0usize;
    let mut skipped_solids = 0usize;
    let mut fallback_identity_solids = 0usize;
    let mut solids_without_shape_representation = 0usize;

    for (data_section, section) in exchange.data.iter().enumerate() {
        let index = step_graph::build_index(&section.entities);
        let document_section = StepSection::new(data_section, &section.entities, &index);
        let ownership = representation_ownership(&section.entities);
        let mut unit_scale_cache = HashMap::<u64, Option<f64>>::new();
        for entity in &section.entities {
            let Some(record) = step_graph::simple_record(entity) else {
                continue;
            };
            if record.name != "MANIFOLD_SOLID_BREP" {
                continue;
            }
            manifold_solids += 1;
            let solid_id = step_graph::entity_id(entity);
            let Some(source) = document_section.entity(solid_id).map(StepEntity::source) else {
                skipped_solids += 1;
                continue;
            };
            debug_assert_eq!(source.data_section, data_section);
            debug_assert_eq!(Some(source.ordinal), index.get(&solid_id).copied());
            let owning = ownership.get(&solid_id).copied();
            let has_shape_representation = owning.is_some();
            solids_without_shape_representation += usize::from(!has_shape_representation);
            let unit_scale_mm = owning.and_then(|(_, context_id)| {
                *unit_scale_cache.entry(context_id).or_insert_with(|| {
                    representation_length_scale_mm(context_id, &document_section)
                })
            });

            let canonical = shape_identity::solid_shape_key(solid_id, &section.entities, &index);
            let (
                geometry_fingerprint,
                identity_kind,
                center_source_units,
                canonical_quarter_turn,
                stats,
            ) = if let Some((key, center, quarter, closure_entities)) = canonical {
                (
                    Some(fingerprint(&key)),
                    IdentityKind::Canonical,
                    Some(center),
                    Some(quarter),
                    BasicSolidStats {
                        closure_entities,
                        vertices: key.vertices,
                        edges: key.edges,
                        oriented_edges: key.oriented_edges,
                        faces: key.faces,
                    },
                )
            } else {
                fallback_identity_solids += 1;
                (
                    None,
                    IdentityKind::FileLocalFallback,
                    None,
                    None,
                    solid_basic_stats(solid_id, &section.entities, &index),
                )
            };

            solids.push(SolidScan {
                data_section: source.data_section,
                solid_id: source.id,
                geometry_fingerprint,
                identity_kind,
                unit_scale_mm,
                center_source_units,
                canonical_quarter_turn,
                closure_entities: stats.closure_entities,
                vertices: stats.vertices,
                edges: stats.edges,
                oriented_edges: stats.oriented_edges,
                faces: stats.faces,
                has_shape_representation,
            });
        }
    }
    solids.sort_by_key(|row| (row.data_section, row.solid_id));
    ScanReport {
        manifold_solids,
        analyzed_solids: solids.len(),
        skipped_solids,
        fallback_identity_solids,
        solids_without_shape_representation,
        solids,
    }
}

#[derive(Debug)]
struct Selection {
    section: usize,
    solid_id: u64,
    fingerprint: String,
}

/// Scan an isolated output with one parse, and require a real product root.
/// A bare shape representation is parseable STEP but not transferable as an
/// independent model by OCCT and other product-oriented STEP importers.
fn scan_isolated_object(bytes: &[u8]) -> Result<ScanReport> {
    let exchange = parse_exchange(bytes)?;
    let report = scan_exchange(&exchange);
    ensure!(
        report.manifold_solids == 1 && report.analyzed_solids == 1,
        "expected exactly one analyzable isolated solid"
    );
    let solid = &report.solids[0];
    let section = exchange
        .data
        .get(solid.data_section)
        .context("missing isolated DATA section")?;
    let ownership = representation_ownership(&section.entities);
    let (representation_id, _) = ownership
        .get(&solid.solid_id)
        .copied()
        .context("isolated solid missing shape representation")?;
    let (_, definition_id) =
        shape_definition_root_for_representation(representation_id, &section.entities)
            .context("isolated solid lacks a transferable product root")?;
    let index = step_graph::build_index(&section.entities);
    ensure!(
        step_entities::named_entity_parameters(
            definition_id,
            "PRODUCT_DEFINITION_SHAPE",
            &section.entities,
            &index
        )
        .is_some(),
        "isolated solid product root has no PRODUCT_DEFINITION_SHAPE"
    );
    Ok(report)
}

fn parse_selection(raw: &str) -> Result<Selection> {
    let mut parts = raw.splitn(3, ':');
    let section = parts
        .next()
        .context("missing DATA section")?
        .parse::<usize>()
        .context("invalid DATA section")?;
    let solid = parts
        .next()
        .context("missing solid id")?
        .parse::<u64>()
        .context("invalid solid id")?;
    let fingerprint = parts.next().context("missing fingerprint")?.to_string();
    if fingerprint.is_empty() || !fingerprint.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        bail!("invalid fingerprint in {raw:?}");
    }
    Ok(Selection {
        section,
        solid_id: solid,
        fingerprint,
    })
}

fn set_shell_faces(entity: &mut EntityInstance, face_ids: &[u64]) -> Result<()> {
    let record = step_graph::simple_record_mut(entity).context("shell is complex")?;
    if record.name != "CLOSED_SHELL" {
        bail!("expected CLOSED_SHELL, got {}", record.name);
    }
    let Parameter::List(params) = &mut record.parameter else {
        bail!("shell parameters are not a list");
    };
    let Some(faces) = params.get_mut(1) else {
        bail!("shell has no face list");
    };
    *faces = Parameter::List(
        face_ids
            .iter()
            .copied()
            .map(step_entities::entity_ref)
            .collect(),
    );
    Ok(())
}

fn set_representation_single_item(entity: &mut EntityInstance, solid_id: u64) -> Result<()> {
    let record = step_graph::simple_record_mut(entity).context("representation is complex")?;
    let Parameter::List(params) = &mut record.parameter else {
        bail!("representation parameters are not a list");
    };
    let Some(items) = params.get_mut(1) else {
        bail!("representation has no item list");
    };
    *items = Parameter::List(vec![step_entities::entity_ref(solid_id)]);
    Ok(())
}

fn extract_one(
    exchange: &Exchange,
    data_section: usize,
    solid_id: u64,
    representation_id: u64,
    context_id: u64,
    index: &HashMap<u64, usize>,
    out_path: &Path,
) -> Result<(String, usize)> {
    let section = exchange
        .data
        .get(data_section)
        .with_context(|| format!("missing DATA section {data_section}"))?;
    let &solid_index = index
        .get(&solid_id)
        .with_context(|| format!("missing solid #{solid_id}"))?;
    let record = step_graph::simple_record(&section.entities[solid_index])
        .with_context(|| format!("solid #{solid_id} is complex"))?;
    if record.name != "MANIFOLD_SOLID_BREP" {
        bail!("#{solid_id} is {}, not MANIFOLD_SOLID_BREP", record.name);
    }

    let shell_id = step_entities::nth_entity_ref(&record.parameter, 1)
        .context("solid has no shell reference")?;
    let face_ids = shape_identity::manifold_solid_face_ids(solid_id, &section.entities, index)
        .context("solid has no semantic face list")?;
    let mut keep = shape_identity::semantic_solid_closure(solid_id, &section.entities, index)
        .context("cannot build semantic solid closure")?;
    keep.extend(step_entities::closure_from(
        context_id,
        &section.entities,
        index,
    ));
    keep.insert(representation_id);

    let mut entities = Vec::with_capacity(keep.len());
    for entity in &section.entities {
        let id = step_graph::entity_id(entity);
        if !keep.contains(&id) {
            continue;
        }
        let mut cloned = entity.clone();
        if id == shell_id {
            set_shell_faces(&mut cloned, &face_ids)?;
        }
        if id == representation_id {
            set_representation_single_item(&mut cloned, solid_id)?;
        }
        entities.push(cloned);
    }

    // Reserve eight fresh, collision-free STEP IDs before adding a self-contained
    // product-definition root for this extracted representation.
    let first_id = section
        .entities
        .iter()
        .map(step_graph::entity_id)
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .context("STEP entity IDs exhausted")?;
    first_id
        .checked_add(7)
        .context("STEP entity IDs exhausted")?;
    entities.extend(isolated_product_root(representation_id, first_id));

    let output = Exchange {
        header: exchange.header.clone(),
        anchor: Vec::new(),
        reference: Vec::new(),
        data: vec![DataSection {
            meta: section.meta.clone(),
            entities,
        }],
        signature: Vec::new(),
    };

    let serialized = crate::write_exchange(&output)?;
    // The isolated graph is already semantics-complete. Run only the normal
    // non-experimental cleanup pass before persisting it: value interning,
    // presentation consolidation and dense entity IDs substantially reduce
    // corpus size without invoking any constructive recovery heuristics.
    let cleaned = crate::clean_bytes(serialized.as_bytes(), &crate::Options::default())?;
    let result = (digest(&cleaned.bytes), cleaned.bytes.len());
    fs::write(out_path, cleaned.bytes).with_context(|| format!("write {}", out_path.display()))?;
    Ok(result)
}

fn extract(bytes: &[u8], out_dir: &Path, selections: &[Selection]) -> Result<Vec<ExtractedBody>> {
    let exchange = parse_exchange(bytes)?;
    fs::create_dir_all(out_dir).with_context(|| format!("create {}", out_dir.display()))?;

    let mut seen = HashSet::new();
    for selection in selections {
        if !seen.insert(&selection.fingerprint) {
            bail!("duplicate output fingerprint {}", selection.fingerprint);
        }
    }

    let mut indexes = Vec::with_capacity(exchange.data.len());
    let mut ownership = Vec::with_capacity(exchange.data.len());
    for section in &exchange.data {
        indexes.push(step_graph::build_index(&section.entities));
        ownership.push(representation_ownership(&section.entities));
    }

    let mut objects = Vec::with_capacity(selections.len());
    for Selection {
        section: section_id,
        solid_id,
        fingerprint,
    } in selections
    {
        let (section_id, solid_id) = (*section_id, *solid_id);
        let index = indexes
            .get(section_id)
            .with_context(|| format!("missing DATA section {section_id}"))?;
        let (representation_id, context_id) = ownership
            .get(section_id)
            .and_then(|owners| owners.get(&solid_id))
            .copied()
            .with_context(|| {
                format!("solid #{solid_id} is not a direct item of a shape representation")
            })?;

        let shard = &fingerprint[..2.min(fingerprint.len())];
        let dir = out_dir.join(shard);
        fs::create_dir_all(&dir).with_context(|| format!("create {}", dir.display()))?;
        let out_path = dir.join(format!("{fingerprint}.step"));
        let (sha256, bytes) = extract_one(
            &exchange,
            section_id,
            solid_id,
            representation_id,
            context_id,
            index,
            &out_path,
        )?;
        objects.push(ExtractedBody {
            fingerprint: fingerprint.clone(),
            sha256,
            bytes,
        });
    }
    Ok(objects)
}

/// Extract selected body fingerprints into standalone STEP files.
/// # Errors
/// Returns an error if selection metadata is invalid or a body cannot be isolated/written.
pub fn extract_selected(
    path: &Path,
    out_dir: &Path,
    selections: &[String],
) -> Result<Vec<ExtractedBody>> {
    let parsed = selections
        .iter()
        .map(|raw| parse_selection(raw))
        .collect::<Result<Vec<_>>>()?;
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    extract(&bytes, out_dir, &parsed)
}

#[cfg(test)]
mod isolated_root_tests {
    use super::*;
    use std::collections::HashSet;

    const BOX: &[u8] = include_bytes!("../validation/fixtures/native_body_corpus_box.step");

    #[test]
    fn every_isolated_body_receives_one_transferable_product_root() -> Result<()> {
        let source = parse_exchange(BOX)?;
        let mut orphaned = source.clone();
        for section in &mut orphaned.data {
            section.entities.retain(|entity| {
                step_graph::simple_record(entity)
                    .is_none_or(|record| record.name != "SHAPE_DEFINITION_REPRESENTATION")
            });
        }
        let orphaned = crate::write_exchange(&orphaned)?.into_bytes();
        // The original body scan succeeds even when OpenCascade has no product
        // to transfer. The extractor must supply the missing root.
        assert!(scan_isolated_object(&orphaned).is_err());
        for input in [BOX, orphaned.as_slice()] {
            let scan = scan_bytes(input)?;
            let solid = scan.solids.first().context("source fixture has no solid")?;
            assert_eq!(scan.analyzed_solids, 1);
            assert!(solid.has_shape_representation);
            let temp = tempfile::tempdir()?;
            let selection = Selection {
                section: solid.data_section,
                solid_id: solid.solid_id,
                fingerprint: "0123456789abcdef".into(),
            };
            let out = temp.path().join("objects");
            let artifacts = extract(input, &out, &[selection])?;
            assert_eq!(artifacts.len(), 1);
            let saved = fs::read(out.join("01/0123456789abcdef.step"))?;
            assert_eq!(artifacts[0].sha256, digest(&saved));
            let exchange = parse_exchange(&saved)?;
            assert_eq!(exchange.data.len(), 1);
            let section = &exchange.data[0];
            let ownership = representation_ownership(&section.entities);
            let output_scan = scan_isolated_object(&saved)?;
            assert_eq!(output_scan.manifold_solids, 1);
            let output_solid = output_scan.solids.first().context("missing output solid")?;
            let (representation_id, _) = ownership
                .get(&output_solid.solid_id)
                .copied()
                .context("isolated representation not found")?;
            let (root, definition) =
                shape_definition_root_for_representation(representation_id, &section.entities)
                    .context("isolated solid lacks a product root")?;
            let index = step_graph::build_index(&section.entities);
            let roots: Vec<_> = section
                .entities
                .iter()
                .filter_map(step_graph::simple_record)
                .filter(|record| record.name == "SHAPE_DEFINITION_REPRESENTATION")
                .collect();
            assert_eq!(roots.len(), 1);
            assert!(index.contains_key(&root));
            let closure = step_entities::closure_from(definition, &section.entities, &index);
            let types: HashSet<_> = closure
                .iter()
                .filter_map(|id| index.get(id))
                .filter_map(|&index| step_graph::simple_record(&section.entities[index]))
                .map(|record| record.name.as_str())
                .collect();
            for kind in [
                "PRODUCT_DEFINITION_SHAPE",
                "PRODUCT_DEFINITION",
                "PRODUCT_DEFINITION_CONTEXT",
                "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE",
                "PRODUCT",
                "PRODUCT_CONTEXT",
                "APPLICATION_CONTEXT",
            ] {
                assert!(types.contains(kind), "missing {kind}");
            }
        }
        Ok(())
    }
}
