//! One parsed STEP exchange with a single typed resolver per DATA section.
//!
//! Import conversion is cached by source-solid ID: topology inspection and
//! editing use the same stable-ID-bearing import, but edits operate on a clone.
//! Neither the original entity graph nor the original imported solid is mutated.

use super::{LoadError, StepSolidImport, Table, convert::StepCompressedTrimmedSolid};
use std::collections::HashMap;
use std::sync::OnceLock;
use step_p21::ast::Exchange;

type CachedImport = OnceLock<Result<StepSolidImport<StepCompressedTrimmedSolid>, String>>;

/// Parsed STEP records, typed entity references and lazily imported B-reps,
/// scoped independently to each DATA section.
#[derive(Debug)]
pub struct ResolvedStepDocument {
    exchange: Exchange,
    tables: Vec<Table>,
    imports: Vec<HashMap<u64, CachedImport>>,
}

impl ResolvedStepDocument {
    /// Resolve every DATA section exactly once.
    #[must_use]
    pub fn from_exchange(exchange: Exchange) -> Self {
        let tables: Vec<Table> = exchange.data.iter().map(Table::from_data_section).collect();
        let imports = tables
            .iter()
            .map(|table| {
                table
                    .manifold_solid_brep
                    .keys()
                    .map(|&id| (id, OnceLock::new()))
                    .collect()
            })
            .collect();
        Self {
            exchange,
            tables,
            imports,
        }
    }

    /// The original STEP AST, including its DATA-section namespaces.
    #[must_use]
    pub const fn exchange(&self) -> &Exchange { &self.exchange }

    /// The typed resolver for one DATA section, if present.
    #[must_use]
    pub fn table(&self, data_section: usize) -> Option<&Table> { self.tables.get(data_section) }

    /// Import a source solid at most once, keeping its retained stable IDs and
    /// load-loss report available to all later read-only or editing operations.
    ///
    /// The import is immutable; callers that apply topology patches must clone
    /// its solid and never edit the cached original.
    ///
    /// # Errors
    /// Returns an error if the source solid is unavailable or its exact trimmed
    /// topology cannot be constructed. Conversion failures are cached as well.
    pub fn import_solid(
        &self,
        data_section: usize,
        solid_id: u64,
    ) -> Result<&StepSolidImport<StepCompressedTrimmedSolid>, LoadError> {
        let table = self
            .tables
            .get(data_section)
            .ok_or_else(|| LoadError::Conversion(format!("missing DATA section {data_section}")))?;
        let solid = table.manifold_solid_brep.get(&solid_id).ok_or_else(|| {
            LoadError::Conversion(format!(
                "missing solid #{solid_id} in DATA section {data_section}"
            ))
        })?;
        let slot = self.imports[data_section].get(&solid_id).ok_or_else(|| {
            LoadError::Conversion(format!(
                "missing import cache for solid #{solid_id} in DATA section {data_section}"
            ))
        })?;
        match slot.get_or_init(|| {
            table
                .to_compressed_trimmed_solid_import(solid)
                .map_err(|error| error.to_string())
        }) {
            Ok(import) => Ok(import),
            Err(detail) => Err(LoadError::Conversion(detail.clone())),
        }
    }

    /// Number of resolved STEP DATA sections.
    #[must_use]
    pub fn len(&self) -> usize { self.tables.len() }

    /// Whether the document contains no DATA sections.
    #[must_use]
    pub fn is_empty(&self) -> bool { self.tables.is_empty() }
}

impl AsRef<Exchange> for ResolvedStepDocument {
    fn as_ref(&self) -> &Exchange { self.exchange() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use step_p21::parser::parse;

    #[test]
    fn section_local_namespaces_are_independent() {
        let text = "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('x'),'1');\nFILE_NAME('a','b',(''),(''),'x','y','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(1.,2.,3.));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(4.,5.,6.));\nENDSEC;\nEND-ISO-10303-21;";
        let document = ResolvedStepDocument::from_exchange(parse(text).unwrap());
        assert_eq!(document.len(), 2);
        assert_eq!(document.exchange().data.len(), 2);
        assert!(document.table(0).unwrap().cartesian_point.contains_key(&1));
        assert!(document.table(1).unwrap().cartesian_point.contains_key(&1));
        assert!(document.table(2).is_none());
    }

    #[test]
    fn one_import_is_reused_without_changing_stable_identity() {
        let source = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../resources/step/occt-cube.step",
        ));
        let document = ResolvedStepDocument::from_exchange(parse(source).unwrap());
        let solid_id = *document
            .table(0)
            .unwrap()
            .manifold_solid_brep
            .keys()
            .next()
            .unwrap();
        let first = document.import_solid(0, solid_id).unwrap();
        let second = document.import_solid(0, solid_id).unwrap();
        assert!(std::ptr::eq(first, second));
        assert_eq!(first.provenance, second.provenance);
        assert!(document.import_solid(1, solid_id).is_err());
    }
}
